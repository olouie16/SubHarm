/*
  ==============================================================================

	This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
SubHarmAudioProcessor::SubHarmAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
	: FFT(fftOrder),
	windowFunc(fftSize, juce::dsp::WindowingFunction<float>::WindowingMethod::hamming),
	AudioProcessor(BusesProperties()
#if ! JucePlugin_IsMidiEffect
#if ! JucePlugin_IsSynth
		.withInput("Input", juce::AudioChannelSet::mono(), true)
#endif
		.withOutput("Output", juce::AudioChannelSet::mono(), true)
#endif
	)
#endif
{
}

SubHarmAudioProcessor::~SubHarmAudioProcessor()
{
}

//==============================================================================
const juce::String SubHarmAudioProcessor::getName() const
{
	return JucePlugin_Name;
}

bool SubHarmAudioProcessor::acceptsMidi() const
{
#if JucePlugin_WantsMidiInput
	return true;
#else
	return false;
#endif
}

bool SubHarmAudioProcessor::producesMidi() const
{
#if JucePlugin_ProducesMidiOutput
	return true;
#else
	return false;
#endif
}

bool SubHarmAudioProcessor::isMidiEffect() const
{
#if JucePlugin_IsMidiEffect
	return true;
#else
	return false;
#endif
}

double SubHarmAudioProcessor::getTailLengthSeconds() const
{
	return 0.0;
}

int SubHarmAudioProcessor::getNumPrograms()
{
	return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
	// so this should be at least 1, even if you're not really implementing programs.
}

int SubHarmAudioProcessor::getCurrentProgram()
{
	return 0;
}

void SubHarmAudioProcessor::setCurrentProgram(int index)
{
}

const juce::String SubHarmAudioProcessor::getProgramName(int index)
{
	return {};
}

void SubHarmAudioProcessor::changeProgramName(int index, const juce::String& newName)
{
}

//==============================================================================
void SubHarmAudioProcessor::prepareToPlay(double sampleRate, int samplesPerBlock)
{
	// Use this method as the place to do any pre-playback
	// initialisation that you need..


	outputRingBuffer.resize(samplesPerBlock + fftSize);
	std::fill(outputRingBuffer.begin(), outputRingBuffer.end(), 0);
	ringBufWriteIndex = 0;
	ringBufReadIndex = samplesPerBlock;

	std::fill(fifo[0], fifo[1] + fftSize, 0);

	subPhase.resize(getNumInputChannels() * sizeof(double));
	std::fill(subPhase.begin(), subPhase.end(), 0);

	m_sampleRate = sampleRate;

	phaseDelta = sampleRate / juce::MathConstants<double>::twoPi;

	subBuffer = juce::AudioBuffer<float>(1, samplesPerBlock);
	subPartialsAmps.resize(nHarmonics, 0);

	//random values for now
	subPartialsPhaseOffsets = { 0.0348901598, 0.9263842084, 0.4589420298, 0.3256211977, 0.8248090600, 0.0418049707, 0.3231477172, 0.2775656167, 0.8058778928, 0.2377176142, 0.4682824886, 0.1786204938, 0.0600388029, 0.1934818492, 0.2626610311, 0.8975879204, 0.5417323311, 0.3042762710, 0.6996761233, 0.5879498383, 0.7859338867, 0.5089617664, 0.4362164675, 0.1965364448, 0.1082312682, 0.5260197481, 0.5935472840, 0.8347254645, 0.5929291021, 0.9293294994 };

	//windowFunc = juce::dsp::WindowingFunction<float>(fftSize, juce::dsp::WindowingFunction<float>::WindowingMethod::hamming);
}

void SubHarmAudioProcessor::releaseResources()
{
	// When playback stops, you can use this as an opportunity to free up any
	// spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool SubHarmAudioProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
#if JucePlugin_IsMidiEffect
	juce::ignoreUnused(layouts);
	return true;
#else
	// This is the place where you check if the layout is supported.
	// In this template code we only support mono or stereo.
	// Some plugin hosts, such as certain GarageBand versions, will only
	// load plugins that support stereo bus layouts.
	if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
		&& layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
		return false;

	// This checks if the input layout matches the output layout
#if ! JucePlugin_IsSynth
	if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
		return false;
#endif

	return true;
#endif
}
#endif

void SubHarmAudioProcessor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midiMessages)
{
	juce::ScopedNoDenormals noDenormals;
	auto totalNumInputChannels = getTotalNumInputChannels();
	auto totalNumOutputChannels = getTotalNumOutputChannels();



	//additive Synthesis Approach
	for (int channel = 0; channel < totalNumOutputChannels; ++channel)
	{
		auto* channelData = buffer.getWritePointer(channel);
		//auto* subBufferData = subBuffer.getWritePointer(0);

		float subSample;

		for (int sample = 0; sample < buffer.getNumSamples(); sample++) {

			channelData[sample] *= 10;//my own mic is to quite;
			pushNextSampleInFifo(channelData[sample], channel);
			if (subActive) {
				subSample = 0;
				for (int partial = 0; partial < subPartialsAmps.size(); partial++) {
					subSample += sin(subPhase[channel] * (2 * partial + 1) + subPartialsPhaseOffsets[partial]) * subPartialsAmps[partial];
				}
				subPhase[channel] += subFreq / phaseDelta;

				channelData[sample] += subSample * volumeOfSlider;
			}
		}


	}

	/*
	* FFT APPROACH
	for (int channel = 0; channel < 1; ++channel)
	{
		auto* channelData = buffer.getWritePointer (channel);
		auto* subBufferData = subBuffer.getWritePointer(0);

		// ..do something to the data...
		for (int sample = 0; sample < buffer.getNumSamples(); sample++) {

			//channelData[sample] *= volumeOfSlider;
			subBufferData[sample] = sin(subPhase[channel]) * 0.1;
			subPhase[channel] += subFreq / phaseDelta;

			//DBG(channelData[sample]);
			//lastSubPhase = subPhase[channel];
			//DBG((subFreq / getSampleRate()) * juce::MathConstants<double>::twoPi);
			//DBG(subPhase);


			pushNextSampleInFifo(channelData[sample], channel, subBufferData[sample]);

		}
		for (int sample = 0; sample < buffer.getNumSamples(); sample++) {
			channelData[sample] = outputRingBuffer[ringBufReadIndex++];
			if (ringBufReadIndex > outputRingBuffer.size()-1)
				ringBufReadIndex = 0;
		}

		if (channel == 1) {
			DBG("#ERROR STEREO CHANNEL!!###");
		}

	}
	*/
}

//==============================================================================
bool SubHarmAudioProcessor::hasEditor() const
{
	return true; // (change this to false if you choose to not supply an editor)
}

juce::AudioProcessorEditor* SubHarmAudioProcessor::createEditor()
{
	return new SubHarmAudioProcessorEditor(*this);
}

//==============================================================================
void SubHarmAudioProcessor::getStateInformation(juce::MemoryBlock& destData)
{
	// You should use this method to store your parameters in the memory block.
	// You could do that either as raw data, or use the XML or ValueTree classes
	// as intermediaries to make it easy to save and load complex data.
}

void SubHarmAudioProcessor::setStateInformation(const void* data, int sizeInBytes)
{
	// You should use this method to restore your parameters from this memory block,
	// whose contents will have been created by the getStateInformation() call.
}

//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
	return new SubHarmAudioProcessor();
}



void SubHarmAudioProcessor::pushNextSampleInFifo(float audioInputSample, int channel, float subOnlySample) {
	if (fifoIndex == fftSize)
	{

		juce::zeromem(fftData, sizeof(fftData));
		memcpy(fftData, fifo[channel], sizeof(fifo[channel]));

		//juce::zeromem(subOnlyFFTData, sizeof(subOnlyFFTData));
		//memcpy(subOnlyFFTData, subFifo, sizeof(subFifo));


		calcFFT(channel);

		fifoIndex = 0;
	}

	fifo[channel][fifoIndex++] = audioInputSample;
}

void SubHarmAudioProcessor::calcFFT(int channel) {

	std::complex<float> freqDomainData[fftSize];


	float mags[fftSize];
	float magsCopy[fftSize];
	float magsPeaksOnly[fftSize];
	float phases[fftSize];
	float intpMags[intpSize];

	int indexToChange;
	float maxMag;


	if (IFFTAproach) {
		std::complex<float> subFreqDomainData[fftSize];
		float subMags[fftSize];
		float subPhases[fftSize];

		for (int i = 0; i < fftSize; i++) {
			complexFFTData[i] = std::complex<float>(fftData[i], 0);
			complexSubFFTData[i] = std::complex<float>(subOnlyFFTData[i], 0);
		}

		FFT.perform(complexFFTData, freqDomainData, false);
		FFT.perform(complexSubFFTData, subFreqDomainData, false);

		for (int i = 0; i < fftSize; i++) {
			mags[i] = std::abs(freqDomainData[i]);
			phases[i] = std::arg(freqDomainData[i]); //in radians
			//subMags[i] = std::abs(subFreqDomainData[i]); not needed
			subPhases[i] = std::arg(subFreqDomainData[i]); //in radians
		}

		//actual calculations in frequency domain missing yet
		//DBG("test");


		//indexToChange = round(subFreq / getSampleRate() * fftSize);
		indexToChange = freq2Index(subFreq, getSampleRate(), fftSize / 2);

		//DBG(indexToChange);

		maxMag = juce::findMaximum(mags, fftSize);
		//DBG(maxMag);

		mags[indexToChange] = 0.2 * fftSize;
		phases[indexToChange] = subPhases[indexToChange];

		mags[indexToChange + 1] = 0.15 * fftSize;
		phases[indexToChange + 1] = subPhases[indexToChange];

		//DBG(phases[indexToChange]-lastPhase);
		//lastPhase = phases[indexToChange];
		//subPhase[channel] += (subFreq / getSampleRate()) * juce::MathConstants<double>::twoPi * fftSize;

		//while (subPhase > juce::MathConstants<float>::twoPi)
		//    subPhase -= juce::MathConstants<float>::twoPi;
		//DBG(subFreq / getSampleRate() * juce::MathConstants<float>::twoPi * fftSize);
		//DBG(subPhase);

		for (int i = 0; i < fftSize; i++) {
			freqDomainData[i] = std::polar(mags[i], phases[i]);
		}

		FFT.perform(freqDomainData, complexFFTData, true); //IFFT


		addToOutputRingBuffer();

	}
	else {//additive Synthesis Approach

		//PitchDetection


		windowFunc.multiplyWithWindowingTable(fftData, fftSize);
		for (int i = 0; i < fftSize; i++) {
			complexFFTData[i] = std::complex<float>(fftData[i], 0);
		}
		FFT.perform(complexFFTData, freqDomainData, false);

		for (int i = 0; i < fftSize; i++) {
			mags[i] = std::abs(freqDomainData[i]);
			phases[i] = std::arg(freqDomainData[i]); //in radians
			magsCopy[i] = mags[i];
		}


		//float estimatedFreq = peakPitchDetection(magsCopy);

		//float estimatedFreq = cepstrumPitchDetection(mags);

		//copyOnlyPeaksFromArray(mags, magsPeaksOnly, fftSize, 2);
		//float estimatedFreq = harmSumPitchDetection(magsPeaksOnly, 1);



		float estimatedFreq = combPitchDetection(fftData);


		DBG(estimatedFreq);


		subFreq = estimatedFreq / 2;


		subActive = estimatedFreq>0;

		//DBG(static_cast<int>(subActive));

		//DBG(maxIndex << ": " << index2Freq(maxIndex, getSampleRate(), fftSize/2) << " / " << estimatedFreq);
		//DBG(mags[maxIndex - 1] << " " << mags[maxIndex] << " " << mags[maxIndex + 1] << " " << estimatedMag);



		//getting mags of harmonics
		std::vector<float> harmAmps = std::vector<float>(nHarmonics, 0);
		std::tuple<double, double> vertexTup;

		for (int harmonic = 0; harmonic < nHarmonics; harmonic++) {
			int harmIndex = freq2Index((harmonic + 1) * estimatedFreq, getSampleRate(), fftSize);

			if (harmIndex > fftSize - 2)
				break;

			vertexTup = calcParabolaVertex(harmIndex - 1, mags[harmIndex - 1], harmIndex, mags[harmIndex], harmIndex + 1, mags[harmIndex + 1]);

			harmAmps[harmonic] = std::get<1>(vertexTup) / fftSize;
		}

		subPartialsAmps[0] = harmAmps[0];
		subPartialsAmps[nHarmonics - 1] = harmAmps[nHarmonics - 1];
		for (int harmonic = 1; harmonic < nHarmonics - 1; harmonic++) {
			subPartialsAmps[harmonic] = (harmAmps[harmonic - 1] + harmAmps[harmonic]) / 2;//maybe log instead of linear middle?
		}



	}
}

void SubHarmAudioProcessor::addToOutputRingBuffer() {

	for (int i = 0; i < fftSize; i++) {
		outputRingBuffer[ringBufWriteIndex++] = complexFFTData[i].real();
		//outputRingBuffer[ringBufWriteIndex++] = fftData[i];
		if (ringBufWriteIndex > outputRingBuffer.size() - 1)
			ringBufWriteIndex = 0;
	}
}



/// <summary>
/// calculates the Vertex of the Parabola going through 3 points
/// </summary>
/// <param name="x1">x of first Point</param>
/// <param name="y1">y of first Point</param>
/// <param name="x2">x of second Point</param>
/// <param name="y2">y of second Point</param>
/// <param name="x3">x of third Point</param>
/// <param name="y3">y of third Point</param>
/// <returns>tuple with (x,y) of Vertex point</returns>
std::tuple<double, double> SubHarmAudioProcessor::calcParabolaVertex(double x1, double y1, double x2, double y2, double x3, double y3) {
	double xv, yv;

	double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
	double A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
	double B = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) / denom;
	double C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

	xv = -B / (2 * A);
	yv = C - B * B / (4 * A);

	return std::make_tuple(xv, yv);
}


/// <summary>
/// return index of maximum value.
/// </summary>
/// <param name="array">array which holds all values</param>
/// <param name="arraySize">size of array</param>
/// <returns>index of maximum value</returns>
int SubHarmAudioProcessor::getIndexOfMaximum(float* array, int arraySize) {
	float maxValue = array[0];
	int maxIndex = 0;
	for (int i = 0; i < arraySize; i++) {
		if (maxValue < array[i]) {
			maxValue = array[i];
			maxIndex = i;
		}
	}
	return maxIndex;
}


int SubHarmAudioProcessor::findNearestIndex(float value, std::vector<float>* vec, int startIndex) {
	//find nearest index to value
	int firstHigherIndex = vec->size() - 1;
	for (int index = startIndex; index < vec->size(); index++) {

		if ((*vec)[index] >= value) {
			firstHigherIndex = index;
			break;
		}
	}
	int nearestIndex = value - (*vec)[firstHigherIndex - 1] < (*vec)[firstHigherIndex] - value && startIndex < firstHigherIndex ? firstHigherIndex - 1 : firstHigherIndex;

	return nearestIndex;
}



/// <summary>
/// copies all relative local peaks with +-peakWidth to left and right to outputArray, while ignoring everything else.
/// Make sure initilizing outputArray before.
/// </summary>
/// <param name="inputArray">array to read peaks from </param>
/// <param name="outputArray">array to copy peaks to</param>
/// <param name="size">size of input and output arrays</param>
/// <param name="peakWidth">number of samples to each side copied with</param>
void SubHarmAudioProcessor::copyOnlyPeaksFromArray(float* inputArray, float* outputArray, int size, int peakWidth) {
	bool rising = inputArray[peakWidth] < inputArray[peakWidth + 1];
	for (int i = peakWidth; i < size - peakWidth; i++) {
		if (rising) {
			if (inputArray[i] > inputArray[i + 1]) {
				//peak detected
				rising = false;
				for (int j = i - peakWidth; j <= i + peakWidth; j++) {
					outputArray[j] = inputArray[j];
				}
			}
		}
		else {
			if (inputArray[i] < inputArray[i + 1]) {
				rising = true;
			}
		}
	}
}


float SubHarmAudioProcessor::cepstrumPitchDetection(float mags[fftSize], float minFreq, float maxFreq) {

	std::complex<float> frequencyDomainData[fftSize];
	std::complex<float> quefrencyDomainData[fftSize];

	for (int i = 0; i < fftSize; i++) {
		//DBG(i << ":" << mags[i] << "-" << log(mags[i]));
		if (mags[i] <= 0) {
			frequencyDomainData[i] = std::complex<float>(-std::numeric_limits<float>::max(), 0);
			DBG("log(0) prevented for: " << mags[i]);
			if (i > 1) {
				DBG("last mag" << mags[i - 1]);
				frequencyDomainData[i] = frequencyDomainData[i - 1] - std::complex<float>(2, 0);
			}
			else {
				DBG("next mag: " << mags[i + 1]);
				frequencyDomainData[i] = std::complex<float>(log(pow(mags[i + 1] / 10, 2)), 0);
			}
		}
		else {
			frequencyDomainData[i] = std::complex<float>(log(pow(mags[i], 2)), 0); //or log10 for base 10 instead of e
		}
	}

	FFT.perform(frequencyDomainData, quefrencyDomainData, false);

	for (int i = 0; i < fftSize; i++) {
		quefrencyDomainData[i] = quefrencyDomainData[i];
	}


	//float dt = 1.0 / (fftSize * (1.0 / getSampleRate()));

	//float df = 1.0 / (fftSize * dt); //

	//int minIndex = ceil((1 / minFreq) / df);
	//int maxIndex = floor((1 / maxFreq) / df);

	int maxIndex = std::min(static_cast<int>(ceil(getSampleRate() / minFreq)), fftSize - 1);
	int minIndex = floor(getSampleRate() / maxFreq);


	auto peakPtr = std::max_element(quefrencyDomainData + minIndex, quefrencyDomainData + maxIndex,
		[](std::complex<float> a, std::complex<float> b) {return std::abs(a) < std::abs(b); } //not sure if real or abs or smt else
	);
	auto minPeakPtr = std::min_element(quefrencyDomainData + minIndex, quefrencyDomainData + maxIndex,
		[](std::complex<float> a, std::complex<float> b) {return std::abs(a) < std::abs(b); } //not sure if real or abs or smt else
	);

	int peakIndex = std::distance(quefrencyDomainData, peakPtr);
	int minPeakIndex = std::distance(quefrencyDomainData, minPeakPtr);

	std::complex<float> peak = *peakPtr;
	std::complex<float> minPeak = *minPeakPtr;
	DBG("max: " << std::real(peak) << "," << std::imag(peak) << ";" << std::abs(peak));
	DBG("min: " << std::real(minPeak) << "," << std::imag(minPeak) << ";" << std::abs(minPeak));

	std::complex<float> average = std::complex<float>(0, 0);
	float averageAbs = 0;
	for (int i = 0; i < fftSize / 2; i++) {
		average += quefrencyDomainData[i];
		averageAbs += abs(quefrencyDomainData[i]);
	}
	average /= fftSize;
	averageAbs /= fftSize;

	DBG("avg: " << std::real(average) << "," << std::imag(average) << "#" << averageAbs);

	std::tuple<double, double> vertexTup = calcParabolaVertex(peakIndex - 1, std::abs(quefrencyDomainData[peakIndex - 1]), peakIndex, std::abs(quefrencyDomainData[peakIndex]), peakIndex + 1, std::abs(quefrencyDomainData[peakIndex + 1]));

	DBG(minIndex << ";" << maxIndex << ";" << std::get<0>(vertexTup));
	return getSampleRate() / std::get<0>(vertexTup);

}




/// <summary>
/// doesn't work properly yet
/// calculates pitch by finding peaks in magnitudes
/// </summary>
/// <param name="mags">array of magnitudes from fft result</param>
/// <returns>pitch</returns>
float SubHarmAudioProcessor::peakPitchDetection(float mags[fftSize]) {

	//copy mags for removing peaks later
	float magsCopy[fftSize];
	std::copy(mags, mags + fftSize, magsCopy);


	//find Peaks in mags
	std::vector<std::tuple<float, float>> peaks; //freq,mag
	float threshhold = 0.01;
	int maxIndex;
	std::tuple<double, double> vertexTup;
	float estimatedMag;
	float estimatedFreq;
	float highestMag;
	int count = 0;
	int maxCount = 50;//to prevent endless loop
	bool searching = true;
	//DBG("_______________");
	while (searching) {
		maxIndex = std::distance(magsCopy, std::max_element(magsCopy, magsCopy + fftSize));

		vertexTup = calcParabolaVertex(maxIndex - 1, magsCopy[maxIndex - 1], maxIndex, magsCopy[maxIndex], maxIndex + 1, magsCopy[maxIndex + 1]);

		estimatedFreq = index2Freq(std::get<0>(vertexTup), getSampleRate(), fftSize / 2);
		estimatedMag = std::get<1>(vertexTup);

		//DBG(estimatedFreq << " - " << estimatedMag);
		//DBG(count);

		peaks.push_back(std::make_tuple(estimatedFreq, estimatedMag));
		if (std::get<1>(peaks[0]) * threshhold > std::get<1>(peaks.back()) || count++ >= maxCount) {
			searching = false;
		}

		//remove Peak
		int leftBorder = 0;
		for (int i = maxIndex - 1; i > 0; i--) {
			if (magsCopy[i] > magsCopy[i + 1]) {
				leftBorder = i;
				break;
			}
		}
		int rightBorder = fftSize - 1;
		for (int i = maxIndex + 1; i <= fftSize - 1; i++) {
			if (magsCopy[i] < magsCopy[i + 1]) {
				rightBorder = i;
				break;
			}
		}
		for (int i = leftBorder; i <= rightBorder; i++) {
			magsCopy[i] = 0;
		}



	}


	//check for harmonic series in Peaks
	/*
	std::sort(peaks.begin(), peaks.end(), [](std::tuple<float, float> a, std::tuple<float, float> b) {return std::get<0>(a) < std::get<0>(b); });
	std::vector<float> peakFreqs(peaks.size());
	for (int i = 0; i < peaks.size(); i++) {
		peakFreqs[i] = std::get<0>(peaks[i]);
	}
	DBG("_________");

	float tolerance = 1 * getSampleRate() / fftSize / 2;

	DBG(tolerance);
	std::string debugStr;

	std::vector<int> foundHarmonics(peakFreqs.size(), 0);
	for (int baseFreqIndex = 0; baseFreqIndex < peakFreqs.size(); baseFreqIndex++) {
		debugStr += std::to_string(peakFreqs[baseFreqIndex]);
		debugStr += "; ";
		int maxHarm = std::min(static_cast<int>(ceil(peakFreqs.back() / peakFreqs[baseFreqIndex])), 10);
		for (int harm = 2; harm < maxHarm; harm++) {
			float calcHarmFreq = harm * peakFreqs[baseFreqIndex];


			int nearestIndex = findNearestIndex(calcHarmFreq, &peakFreqs, baseFreqIndex+1);
			DBG("nearest: " << calcHarmFreq << "|" << peakFreqs[nearestIndex]);
			//check if in tolerance
			//DBG((peakFreqs[baseFreqIndex] + tolerance)* harm << ">" << peakFreqs[nearestIndex]);
			//DBG((peakFreqs[baseFreqIndex] - tolerance)* harm << "<" << peakFreqs[nearestIndex]);
			if ((peakFreqs[baseFreqIndex] + tolerance) * harm > peakFreqs[nearestIndex] && (peakFreqs[baseFreqIndex] - tolerance) * harm < peakFreqs[nearestIndex]){
				//in tolerance
				foundHarmonics[baseFreqIndex]++;
				DBG("harms: " << peakFreqs[baseFreqIndex] << ": " << harm << " - " << peakFreqs[nearestIndex]);
			}
			else {
				// -- as counter for very low freqs having almost everything in tolerance, but not every harmonic is present
				foundHarmonics[baseFreqIndex]--;
			}
		}
	}
	maxIndex = std::distance(foundHarmonics.begin(), std::max_element(foundHarmonics.begin(), foundHarmonics.end()));
	estimatedFreq = peakFreqs[maxIndex];
	DBG(debugStr);
	DBG(estimatedFreq << ":" << foundHarmonics[maxIndex]);
	*/

	return std::get<0>(peaks[0]);

}


/// <summary>
/// calculates pitch by summing the harmonics. If no or too quite harmonic Sound, returns -1.
/// </summary>
/// <param name="mags">magnitudes Spectrum eg by abs(fftResult) </param>
/// <param name="nHarmonics">number of Harmonics taken into account</param>
/// <returns>frequency; or if no harmonic sound was found -1</returns>
float SubHarmAudioProcessor::harmSumPitchDetection(float mags[fftSize], int nHarmonics) {

	float sumMags[fftSize];
	float avg = 0;
	//DBG("________________________");
	for (int sample = 0; sample < fftSize; sample++) {
		sumMags[sample] = 0;
		for (int harmonic = 1; harmonic <= std::min(nHarmonics, fftSize / (sample + 1)); harmonic++) {
			//DBG(sample << "*" << harmonic << "=" << sample * harmonic);
			sumMags[sample] += mags[sample * harmonic] / (100 + harmonic) * 100;
		}
		avg += sumMags[sample];
	}
	avg /= fftSize;
	//DBG("________________________");

	int peakIndex = getIndexOfMaximum(sumMags, fftSize);
	std::tuple<double, double> vertexTup = calcParabolaVertex(peakIndex - 1, sumMags[peakIndex - 1], peakIndex, sumMags[peakIndex], peakIndex + 1, sumMags[peakIndex + 1]);
	float freq = index2Freq(std::get<0>(vertexTup), getSampleRate(), fftSize / 2);

	//DBG("avg: " << avg);
	//DBG("max: " << std::get<1>(vertexTup));

	//return -1 if no or to quite harmonic Sound
	if (std::get<1>(vertexTup) / avg > 100) {
		return freq;
	}
	else {
		return -1;
	}
}



float SubHarmAudioProcessor::combPitchDetection(float samples[fftSize], float minFreq, float maxFreq) {
	
	/*
	combfilter by adding delayed signal. 
	delay=wavelenght/2 of first combfreq
	

	*/



	int maxDelay = getSampleRate() / minFreq / 2; //in samples
	int minDelay = getSampleRate() / maxFreq / 2;

	std::vector<float> rms = std::vector<float>(maxDelay - minDelay+1);
	float filteredSamples[fftSize];
	float rmsAvg = 0;

	for (int delay = minDelay; delay <= maxDelay; delay++) {

		rms[delay - minDelay] = 0;

		//calc buffer with delayed signal
		for (int sample = 0; sample < fftSize; sample++) {
			//could be optimized
			//if (sample <= delay) {
			//	filteredSamples[sample] = samples[sample];
			//}
			//else {
			//	filteredSamples[sample] = samples[sample] + samples[sample - delay];
			//}

			//rms[delay - minDelay] += pow(filteredSamples[sample], 2);

			if (sample <= delay) {
				rms[delay - minDelay] += pow(samples[sample], 2);
			}
			else {
				rms[delay - minDelay] += pow(samples[sample] + samples[sample - delay], 2);
			}

		}

		rms[delay - minDelay] = sqrt(rms[delay - minDelay] / fftSize);
		rmsAvg += rms[delay - minDelay];

	}

	rms.empty() ? 0 : rmsAvg /= rms.size();

	int minIndex = std::distance(std::begin(rms), std::min_element(std::begin(rms), std::end(rms)));
	int maxIndex = std::distance(std::begin(rms), std::max_element(std::begin(rms), std::end(rms)));

	DBG("avg: " << rmsAvg);
	DBG("min: " << rms[minIndex]);
	DBG("max: " << rms[maxIndex]);
	float freq = getSampleRate() / (2 * (minIndex + minDelay));
	DBG("freq: " << freq);

	return freq;
}
