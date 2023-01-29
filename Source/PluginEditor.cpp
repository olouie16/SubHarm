/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
SubHarmAudioProcessorEditor::SubHarmAudioProcessorEditor (SubHarmAudioProcessor& p)
    : AudioProcessorEditor (&p), audioProcessor (p)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    setSize (400, 300);

    volumeSlider.setSliderStyle(juce::Slider::LinearBarVertical);
    // these define the parameters of our slider object
    volumeSlider.setSliderStyle(juce::Slider::LinearBarVertical);
    volumeSlider.setRange(0.0, 5.0);
    volumeSlider.setTextBoxStyle(juce::Slider::NoTextBox, false, 90, 0);
    volumeSlider.setPopupDisplayEnabled(true, false, this);
    volumeSlider.setTextValueSuffix(" Volume");
    volumeSlider.setValue(1.0);

    // this function adds the slider to the editor
    addAndMakeVisible(&volumeSlider);

    // add the listener to the slider
    volumeSlider.addListener(this);

}

SubHarmAudioProcessorEditor::~SubHarmAudioProcessorEditor()
{
}

//==============================================================================
void SubHarmAudioProcessorEditor::paint (juce::Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));

    g.setColour (juce::Colours::white);
    g.setFont (15.0f);
    g.drawFittedText ("Hello World!2", getLocalBounds(), juce::Justification::centred, 1);
}

void SubHarmAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..

    volumeSlider.setBounds(40, 30, 20, getHeight() - 60);

}
void SubHarmAudioProcessorEditor::sliderValueChanged(juce::Slider* slider)
{
    audioProcessor.volumeOfSlider = volumeSlider.getValue();
}