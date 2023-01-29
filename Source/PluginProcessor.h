/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>

//==============================================================================
/**
*/
class SubHarmAudioProcessor  : public juce::AudioProcessor
                            #if JucePlugin_Enable_ARA
                             , public juce::AudioProcessorARAExtension
                            #endif
{
public:
    //==============================================================================
    SubHarmAudioProcessor();
    ~SubHarmAudioProcessor() override;

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const juce::String getProgramName (int index) override;
    void changeProgramName (int index, const juce::String& newName) override;

    //==============================================================================
    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

    float volumeOfSlider;
    
    const static int fftOrder = 12;
    const static int fftSize = 1 << fftOrder;

    const static int intpRatio = 5;
    const static int intpSize = fftSize * intpRatio;


private:
    //==============================================================================
    
    bool IFFTAproach = false; //true = additive Synth Approach
    
    juce::dsp::FFT FFT;
    juce::dsp::WindowingFunction<float> windowFunc;
    float fifo[2][fftSize];
    float fftData[fftSize];
    float subFifo[fftSize];
    float subOnlyFFTData[fftSize];
    int fifoIndex = 0;
    bool nextFFTBlockReady = false;
    std::complex<float> complexFFTData[fftSize];
    std::complex<float> complexSubFFTData[fftSize];

    std::vector<float> outputRingBuffer;
    int ringBufReadIndex;
    int ringBufWriteIndex;

    //float subFreq = 82.4069; //E2
    float subFreq = 121;
    std::vector<double> subPhase;

    int m_sampleRate;
    float phaseDelta;
    float lastPhase;
    juce::AudioBuffer<float> subBuffer;
    std::vector<float> subPartialsAmps;
    std::vector<float> subPartialsPhaseOffsets;
    int nHarmonics = 30;
    bool subActive = true;

    void pushNextSampleInFifo(float sample, int channel, float subOnlySample=0);
    void calcFFT(int channel);
    void addToOutputRingBuffer();
    float peakPitchDetection(float mags[fftSize]);
    float cepstrumPitchDetection(float mags[fftSize], float minFreq=50, float maxFreq=200);
    float harmSumPitchDetection(float mags[fftSize], int nHarmonics);
    float combPitchDetection(float samples[fftSize], float minFreq=80, float maxFreq=300);

    double index2Freq(double i, double samplerate, int fftSize) {
        return (double)i * (samplerate / fftSize / 2.0);
    }

    int freq2Index(double freq, double samplerate, int fftSize) {
        return (int)(freq / (samplerate / fftSize / 2.0));
    }

 
    std::tuple<double, double> calcParabolaVertex(double x1, double y1, double x2, double y2, double x3, double y3);
    int getIndexOfMaximum(float* array, int arraySize);
    int findNearestIndex(float value, std::vector<float>* vec, int startIndex=0);
    void copyOnlyPeaksFromArray(float* inputArray, float* outputArray, int size, int peakWidth=2);

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (SubHarmAudioProcessor)
};
