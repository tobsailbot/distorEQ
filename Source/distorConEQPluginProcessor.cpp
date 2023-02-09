#include "../JuceLibraryCode/JuceHeader.h"

#include "distorConEQ.h"
#include "distorConEQ_types.h"

#if JUCE_VERSION >= 0x050400
using Parameter = AudioProcessorValueTreeState::Parameter;
#endif

struct onParamChangeListener : AudioProcessorValueTreeState::Listener
{
    onParamChangeListener(distorConEQStackData* sd)
    : SD(sd)
    {
    }
    
    void parameterChanged (const String& parameterID, float newValue) override
    {
	    (void)parameterID;
	    int idx = -1;
        if (parameterID == "gain") {
            idx = 0;
        } else if (parameterID == "fEQ") {
            idx = 1;
        } else if (parameterID == "gainEQ") {
            idx = 2;
        } else if (parameterID == "level") {
            idx = 3;
        }
		onParamChangeCImpl(SD, idx, static_cast<double>(newValue));
    }

    distorConEQStackData *SD;
};

//==============================================================================
class distorConEQAudioProcessor  : public AudioProcessor
{
    //==============================================================================
#if JUCE_VERSION >= 0x050400

public:
    distorConEQAudioProcessor()
        : paramListener(&mStackData),
          parameters(*this, nullptr, "distorConEQ", {
                std::make_unique<Parameter>("gain", "Drive", "",
                    NormalisableRange<float>(1.f,10.f), 1.f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("fEQ", "Frequency", "",
                    NormalisableRange<float>(80.f,16000.f,[](float min, float max, float norm) {return min*powf(max/min,norm);}, [](float min, float max, float val) {return logf(val/min)/logf(max/min);}), 1000.f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("gainEQ", "Frequency gain", "",
                    NormalisableRange<float>(-12.f,12.f), 0.f, [](float val) {return String(val, 3);}, nullptr),
                std::make_unique<Parameter>("level", "Master level", "",
                    NormalisableRange<float>(0.f,2.f), 1.f, [](float val) {return String(val, 3);}, nullptr) })

    {
        mStackData.pd = &mPersistentData;
        
        distorConEQ_initialize(&mStackData);

        createPluginInstance(&mStackData, reinterpret_cast<unsigned long long>(this));

        parameters.addParameterListener("gain", &paramListener);
        parameters.addParameterListener("fEQ", &paramListener);
        parameters.addParameterListener("gainEQ", &paramListener);
        parameters.addParameterListener("level", &paramListener);

    }
    //==============================================================================
#else // For JUCE prior to 5.4.0
public:
    distorConEQAudioProcessor()
    :   paramListener(&mStackData), parameters (*this, nullptr)
    {
        mStackData.pd = &mPersistentData;
        
        distorConEQ_initialize(&mStackData);

        createPluginInstance(&mStackData, reinterpret_cast<unsigned long long>(this));

        //
        // Parameter property gain
        //
        parameters.createAndAddParameter ("gain", "Drive", "",
            NormalisableRange<float>(1.f, 10.f), 1.f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("gain", &paramListener);

        //
        // Parameter property fEQ
        //
        parameters.createAndAddParameter ("fEQ", "Frequency", "",
            NormalisableRange<float>(80.f, 16000.f, 
                [](float min, float max, float norm) {return min * powf(max/min, norm);},
                [](float min, float max, float val) {return logf(val/min)/logf(max/min);} ),
            1000.f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("fEQ", &paramListener);

        //
        // Parameter property gainEQ
        //
        parameters.createAndAddParameter ("gainEQ", "Frequency gain", "",
            NormalisableRange<float>(-12.f, 12.f), 0.f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("gainEQ", &paramListener);

        //
        // Parameter property level
        //
        parameters.createAndAddParameter ("level", "Master level", "",
            NormalisableRange<float>(0.f, 2.f), 1.f,
            [](float val) {return String(val, 3);},
            nullptr);
        parameters.addParameterListener("level", &paramListener);

        parameters.state = ValueTree(Identifier("distorConEQ"));
    }
#endif

    //==============================================================================
    ~distorConEQAudioProcessor()
    {
        distorConEQ_terminate(&mStackData);
    }
    
    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override
    {
        (void)samplesPerBlock;
        resetCImpl(&mStackData, sampleRate);
        setLatencySamples(getLatencyInSamplesCImpl(&mStackData));
    }

    void releaseResources() override                { }
    
    
    void processBlock (AudioBuffer<double>& buffer, MidiBuffer& midiMessages) override
    {
        (void)midiMessages;
        ScopedNoDenormals noDenormals;
        const double** inputs = buffer.getArrayOfReadPointers();
        double** outputs = buffer.getArrayOfWritePointers();
        int nSamples = buffer.getNumSamples();
        distorConEQStackData *SD = &mStackData;

        int osz0_;
        int osz1_;
        if (nSamples <= MAX_SAMPLES_PER_FRAME) {
            /* Fast path for common frame sizes. */
            int isz0_ = nSamples;
            int isz1_ = nSamples;
            processEntryPoint(SD, (double)nSamples,
                    inputs[0], &isz0_,
                    inputs[1], &isz1_,
                    outputs[0], &osz0_,
                    outputs[1], &osz1_);
        } else {
            /* Fallback for unusually large frames. */
            int isz0_ = MAX_SAMPLES_PER_FRAME;
            int isz1_ = MAX_SAMPLES_PER_FRAME;
            int n = MAX_SAMPLES_PER_FRAME;
            for (int i_ = 0; i_ < nSamples; i_ += MAX_SAMPLES_PER_FRAME) {
                if (i_ + MAX_SAMPLES_PER_FRAME > nSamples) {
                    n = nSamples - i_;
                    isz0_ = nSamples - i_;
                    isz1_ = nSamples - i_;
                }
                processEntryPoint(SD, (double)n,
                        inputs[0]+i_, &isz0_,
                        inputs[1]+i_, &isz1_,
                        outputs[0]+i_, &osz0_,
                        outputs[1]+i_, &osz1_);
            }
        }

    }
    
    void processBlock (AudioBuffer<float>& buffer,  MidiBuffer& midiMessages) override
    {
        (void)midiMessages;
        AudioBuffer<double> doubleBuffer;
        doubleBuffer.makeCopyOf(buffer);
        processBlock(doubleBuffer, midiMessages);
        buffer.makeCopyOf(doubleBuffer);
    }
    
    //==============================================================================
    bool hasEditor() const override                 { return true; }
    AudioProcessorEditor* createEditor() override;
    
    //==============================================================================
    const String getName() const override           { return JucePlugin_Name; }

    bool acceptsMidi() const override               { return false; }
    bool producesMidi() const override              { return false; }
    bool isMidiEffect () const override             { return false; }
    double getTailLengthSeconds() const override    { return 0.0;   }

    //==============================================================================
    // NB: some hosts don't cope very well if you tell them there are 0 programs,
    // so this should be at least 1, even if you're not really implementing programs.
    int getNumPrograms() override                       { return 1;  }
    int getCurrentProgram() override                    { return 0;  }
    void setCurrentProgram (int index) override         { (void) index; }
    const String getProgramName (int index) override    { (void) index; return {}; }
    void changeProgramName (int index, const String& newName) override  { (void) index; (void) newName; }
    
    //==============================================================================
    void getStateInformation (MemoryBlock& destData) override
    {
        auto xml (parameters.state.createXml());
        copyXmlToBinary (*xml, destData);
    }
    
    void setStateInformation (const void* data, int sizeInBytes) override
    {
        auto xmlState (getXmlFromBinary (data, sizeInBytes));
        if (xmlState != nullptr)
            if (xmlState->hasTagName (parameters.state.getType()))
                parameters.state = ValueTree::fromXml (*xmlState);
    }
    
    bool supportsDoublePrecisionProcessing() const override  { return true; }
    
private:
    //==============================================================================
    static const int MAX_SAMPLES_PER_FRAME = 4096;

    distorConEQStackData mStackData;
    distorConEQPersistentData mPersistentData;
    onParamChangeListener paramListener;
    
    //==============================================================================
    AudioProcessorValueTreeState parameters;
 
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (distorConEQAudioProcessor)
};

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new distorConEQAudioProcessor();
}

#include "distorConEQPluginEditor.h"

AudioProcessorEditor* distorConEQAudioProcessor::createEditor()
{
    return new distorConEQAudioProcessorEditor(*this, parameters);
}

