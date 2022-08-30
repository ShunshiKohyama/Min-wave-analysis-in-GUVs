// Min_wave_analysis_in_GUVs
// ImageJ macro for detecting GUVs, making kymographs, and anlysing the wave patterns and periods based on fluorescent signals of Min waves inside GUVs
// Validated only for Fiji v1.53f or later
// Input: 2D or 3D-stacked time-lapse image files with fluorescently labeled (inside or membrane) GUVs
// Resulting data sets are saved as new Tiff image files (for kymographs and time-lapse images), a CSV file (for pattern and period analysis), or a Zip file (for ROIs)
// For detailed information: Kohyama et al, Nature Communications 2022
//
// Shunshi Kohyama (Max Planck Institute of Biochemistry), 2022
// Leonard Fröhlich (Max Planck Institute of Biochemistry), 2022

// Set image information and target channel, slice, threshold method, and threshold values
methods = getList("threshold.methods");
getVoxelSize(width, height, depth, unit);
interval =  Stack.getFrameInterval();
Dialog.create("Min wave analysis");
Dialog.addString("Title", getTitle());
if (unit != "µm") {
	Dialog.addMessage("                          Please check!!", 12, "red");
}
Dialog.addNumber("Set scale (µm in a pixel)", width);
if (interval == 0) {
	Dialog.addMessage("                          Please check!!", 12, "red");
}
Dialog.addNumber("Set interval (sec)", interval);
Dialog.addNumber("Set channel (color)", 1);
Dialog.addNumber("Set slice (z position)", 1);
Dialog.addNumber("Set frame (time point)", 1);
Dialog.addChoice("Threshold method", methods, "Huang");
Dialog.addNumber("Minimal diameter for detection (µm)", 10);
Dialog.addNumber("Maximal diameter for detection (µm)", 50);
Dialog.addCheckbox("fast mode", false);
Dialog.show();
title = Dialog.getString();
scale = Dialog.getNumber();
interval = Dialog.getNumber();
channel = Dialog.getNumber();
slice = Dialog.getNumber();
frame = Dialog.getNumber();
type = Dialog.getChoice();
min = Dialog.getNumber();
max = Dialog.getNumber();
fast = Dialog.getCheckbox();
minArea = PI * pow((min/2), 2);
maxArea =  PI * pow((max/2), 2);
if (fast == 1) {
	if (isOpen("ROI Manager")) {
	selectWindow("ROI Manager");
	run("Close");
	}
} else {
	run("ROI Manager...");
}
run("Set Scale...", "distance=1 known=scale pixel=1 unit=µm");

setBatchMode(true);
// Set target channel, slice and then duplicate, apply thresholding, and add ROIs to the ROI manager
run("Duplicate...", "duplicate channels=channel slices=slice");
Stack.setFrame(frame);
run("Duplicate...", " ");
setAutoThreshold(type+" dark");
run("Convert to Mask");
run("Analyze Particles...", "size=minArea-maxArea circularity=1.00 clear include add");
close();
roiManager("Show All with labels");

// Manually delete or add ROIs, the ROIs should be slightly larger circle than the target GUVs peripheral
setBatchMode("show");
beep();
waitForUser("Add or delete ROIs");
Dialog.create("Circle detection");
Dialog.addNumber("Searching range (um)", 5);
Dialog.addCheckbox("Skip detection?", false);
Dialog.show();
range = Dialog.getNumber();
skip = Dialog.getCheckbox();

// Determine the GUVs by scanning the highest weighted-intensity profile of the peripheral
setBatchMode("hide");
N = roiManager("count");
pixels = Math.ceil(range/scale);
for (i=0; i<N; i++) {
	roiManager("select", 0);
	getSelectionBounds(x, y, w, h); // Get the boundary (starting position and size) of the selected ROI
	d = minOf(w, h); // Set the d (miximal diameter of the circle for searching) as smaller length of the either x or y axis
	if (pixels > d) {
		pixels = d;
	}
	roiManager("Delete");
	
	I = 0;
	if (skip == false) {
		for (j=0; j<pixels; j++) {
			makeOval(x, y, d-j, d-j); // Make the circle at top-left position inside selected ROI, the size of the ROI decreases 
			run("Area to Line"); // convert the circle to circular line (as it is for calculation of the peripheral intensity)
			for (k=0; k<=h-d+j; k++) {
				for (l=0; l<=w-d+j; l++) {
					Roi.move(x+l, y+k); // move the circular line to the every position inside the boundary
					getRawStatistics(nPixels, mean, min, max, std, histogram); // get the statistics of the circular line
					if (I < mean * nPixels) { // calculate the weghted-intensity (total intensity) of the peripheral and rewrite the ROI with the highest value
						I = mean * nPixels;
						roiManager("add");
					}
				}
			}
			while (roiManager("count") > N) { // Keep the ROI with the highest intensity and delete the rest
				roiManager("Select", N-1);
				roiManager("Delete");
			}
		}
	} else { // Skip the detection steps
		makeOval(x, y, w, h);
		run("Area to Line");
		roiManager("add");
	}
}
Roi.remove;

// Manually add or delete ROIs (here newly added ROIs must be circular lines but not circles by running the "Area to Line" command)
setBatchMode("show");
Overlay.hide;
roiManager("show all with labels");
beep();
waitForUser("Add or delete ROIs");

// Choose the actions and set searching range for floating calibration
Dialog.create("Analyze setting");
labels = newArray("Wave analysis", "Manual pattern detection", "Show results", "Save results (csv)", "Floating calibration", "Save calibrated stack", "Save kymographs", "Save ROIs");
defaults = newArray(true, false, true, false, false, false, false, false);
answer = newArray(8);
Dialog.addCheckboxGroup(4, 2, labels, defaults);
Dialog.addNumber("Searching range (floating calibration) (µm)", 2);
Dialog.show();
for (i=0; i<8; i++) {
	answer[i] = Dialog.getCheckbox();
}
range = Dialog.getNumber();
n = nSlices;
N = roiManager("count");

if (answer[3] == 1 || answer[5] == 1 || answer[6] == 1 || answer[7] == 1) {
	dir = getDirectory("Choose a Directory");
}
	
if (answer[4] == 1) { // Make arrays for floating calibration
	X = newArray(n);
	Y = newArray(n);
	name = getTitle();
	Overlay.hide;
}

if (answer[0] == 1) {
	smoothingArray = newArray(1, 2, 4, 5, 8, 10, 20, 40);
	Dialog.create("Oscillation analysis");
	Dialog.addRadioButtonGroup("Smoothing parameter", smoothingArray, 2, 4, "10.0");
	Dialog.addNumber("Min. period (s)", 30);
	Dialog.addNumber("Max. period (s)", 180);
	Dialog.addNumber("Min. Angle (°) for traveling", 5);
	Dialog.addNumber("Min. period-ratio for pulsing", 0.9);
	Dialog.addNumber("Max. period-ratio for pulsing", 1.1);
	Dialog.addNumber("Min. period-ratio for pole-to-pole", 1.8);
	Dialog.addNumber("Max. period-ratio for pole-to-pole", 2.2);
	Dialog.show();
	smoothing = Dialog.getRadioButton();
	minPeriod = Dialog.getNumber();
	maxPeriod = Dialog.getNumber();
	minAngle = Dialog.getNumber();
	minRatioPuls = Dialog.getNumber();
	maxRatioPuls = Dialog.getNumber();
	minRatioP2P = Dialog.getNumber();
	maxRatioP2P = Dialog.getNumber();
}

setBatchMode("hide");
run("Clear Results");
if (answer[7] == 1) {
	roiManager("save", dir+title+"_ROI.zip") // Save the set of ROIs
}

for (i=1; i<=N; i++) {
	roiManager("select", i-1);

	// Get the ID and diameter of the selected ROI
	if (answer[0] == 1) { 
		setResult("ID", i-1, title+"_"+i);
		setResult("diameter (µm)", i-1, getValue("Height"));
	}

	// Floating calibration
	if (answer[4] == 1) {
		pixels = Math.ceil(range/scale);
		getSelectionBounds(x, y, w, h); // Get the position of the first (t=0) image
		X[0] = x; // Fill the first position of x (same for y)
		Y[0] = y;

		for (j=2; j<=n; j++) { // Detect the position of the cell from j=2 (t=1) to the finall image of the sequence
			setSlice(j);
			I = 0;
			for (k=-pixels; k<=pixels; k++) {
				for (l=-pixels; l<=pixels; l++) {
					Roi.move(X[j-2]+l, Y[j-2]+k); // move the ROI (which was already choosen as a circular line) around +-S range from the position of the previous image
					getRawStatistics(nPixels, mean, min, max, std, histogram);
					if (I < mean) { // Get the mean intensity of the peripheral and rewrite with the highest value
						I = mean;
						X[j-1] = X[j-2]+l; // and record the current position of x by adding the difference from previous position (same for y)
						Y[j-1] = Y[j-2]+k;
					}
				}
			}
		}
		Overlay.hide;

		// Make a new floating-calibrated time-lapse stack
		for (j=1; j<=n; j++) {
			selectWindow(name); // Select the original image
			setSlice(j);
			makeRectangle(X[j-1]-pixels, Y[j-1]-pixels, w+2*pixels, h+2*pixels); // Recall the calibrated position of the GUV at each time point and make a square ROI to center the GUV
			run("Duplicate...", " "); // Duplicate the ROI, resultant image shows the GUV at the center of the image
		}
		selectWindow(name);
		Roi.remove;
		run("Images to Stack", "name=Stack title=[]"); // Make a new time-lapse stack by all duplicated images
		
		if (answer[5] == 1) {
			saveAs("Tiff", dir+title+"_"+i+".tif");
		}
		if (answer[0] == 1 || answer[6] == 1) {
			makeOval(pixels, pixels, w, h); // In the case of making kymograph, make a circle along with the peripheral and then make a circular line
			run("Area to Line");
		}
	}

	// Make kymographs
	if (answer[0] == 1 || answer[6] == 1) {
		run("Straighten...", "title=straight line=5 process"); // Convert the circular line into the straight line with 5 pixel thickness along with the line
		w = getValue("Width"); // Get the length of the straightened line in pixel
		run("Size...", "width=w height=1 depth=n average interpolation=Bilinear"); // Convert the image size into 1 pixel height (with the same width)
		run("Make Montage...", "columns=1 rows=n scale=1"); // Make a kymograph by periperal intensity (width) vs time-points (height)
		if (answer[6] == 1) {
			saveAs("Tiff", dir+title+"_"+i+"_kymo.tif");
			if (answer[0] == 0) {
				close(); // For kymograph
				close(); // For straight line
			}
		}
	}
	
	// Wave analysis
	if (answer[0] == 1) {
		// Manual pattern detection
		if (answer[1] == 1) {
			setBatchMode("show");
			Dialog.create("Pattern detection");
			patterns = newArray("Pulsing", "Pole-to-pole", "Traveling", "Undefined");
			Dialog.addRadioButtonGroup("patterns", patterns, 1, 4, "Pulsing");
			Dialog.show();
			pattern = Dialog.getRadioButton();
			setBatchMode("hide");
			setResult("Pattern", i-1, pattern);
		}
		
		// Directionality detection
		run("Duplicate...", "title=gx"); // Duplicate the kymograph for 5x5 sobel filtering in x and y derivative
		run("Duplicate...", "title=gy");
		selectWindow("gx");
		run("Convolve...", "text1=[2 2 4 2 2 \n1 1 2 1 1\n0 0 0 0 0\n-1 -1 -2 -1 -1\n-2 -2 -4 -2 -2\n] normalize");
		selectWindow("gy");
		run("Convolve...", "text1=[2 1 0 -1 -2\n2 1 0 -1 -2\n4 2 0 -2 -4\n2 1 0 -1 -2\n2 1 0 -1 -2\n] normalize");
		imageCalculator("Divide create", "gy","gx"); // Generating a new image by calculating gy/gx
		close("gx");
		close("gy");
		selectWindow("Result of gy");
		w = getWidth();
		h = getHeight();
		UArray = newArray(w*h);
		for (y=0; y<h; y++) {
			for (x=0; x<w; x++) {
				UArray[x+(w*y)]= 180 * atan(getPixel(x, y)) / PI; // Obtain the directionality of each pixel by calculating arctan(dx/dy)
			}
		}

		Plot.create("Histogram", "Angle [°]", "Frequency"); //put each pixel value into a histogram
		Plot.setLimits(-90, 90, 0, w*h);
		Plot.addHistogram(UArray, 1, 0);
		Plot.show();
		IJ.renameResults("Results","toSave");
		selectWindow("Histogram");
		Plot.showValues();
		angle = newArray(nResults);
		frequency = newArray(nResults);
		for (rowIndex = 0; rowIndex < nResults; rowIndex++) {
			angle[rowIndex] = getResult("X", rowIndex);
			frequency[rowIndex] = getResult("Y", rowIndex);
		}
		close("Results");
		close("Histogram");
		close("Result of gy");
		Fit.doFit("Gaussian", angle, frequency); // Fit gausian curve to get the most prominent direction
		mainDirection = Fit.p(2);
		run("Clear Results");
		IJ.renameResults("toSave","Results");
		setResult("Main Direction [°]", i-1, mainDirection);
		
		// Period detection
		h = getHeight();
		w = getWidth();
		int = 0;
		pos = 0;
		makeRectangle(0, 0, 5, h); // Choose the most prominent position of oscillation in the kymograph and duplicate
		for (j = 0; j <= w-5; j++) {
			Roi.move(j, 0);
			temp = getValue("Mean");
			if (temp > int) {
				int = temp;
				pos = j;
			}
		}
		makeRectangle(pos, 0, 5, h);
		run("Duplicate...", " ");
		height = h * smoothing; // Interpolate the pixels along with time points by using smoothing factor to increase the data points for fitting
		run("Size...", "width=1 height=height depth=1 average interpolation=Bilinear");
		makeLine(0, 0, 0, height);
		profile = getProfile();
		close(); // For 1d line (from selected position)
		
		timepoints = Array.getSequence(lengthOf(profile)); // Fit Sin wave function to obtain all parameters
		for (j=0; j<lengthOf(profile); j++) {
			timepoints[j] = j * interval / smoothing;
		}
		Fit.doFit("y = a * sin(b * x + c) + d", timepoints, profile);
		R = Fit.rSquared;
		periodRaw = Fit.p(1);
		repeat = Math.floor(periodRaw/PI);
		periodAmp = periodRaw - PI * repeat;
		if (periodAmp > PI/2) {
			periodAmp = PI - periodAmp;
		}
		period = 2 * PI / periodAmp;
		setResult("Period (sec)", i-1, period);
		setResult("Fitting (R^2)", i-1, R);
		Roi.remove;
		
		// Repeat the same sequence but for the entire width of the kymograph
		run("Size...", "width=1 height=height depth=1 average interpolation=Bilinear");
		makeLine(0, 0, 0, height);
		profileCrop = getProfile();
		close(); // For 1d line (from original kymograph)
		Fit.doFit("y = a * sin(b * x + c) + d", timepoints, profileCrop);
		RCrop = Fit.rSquared;
		periodRawCrop = Fit.p(1);
		repeatCrop = Math.floor(periodRawCrop/PI);
		periodAmpCrop = periodRawCrop - PI * repeatCrop;
		if (periodAmpCrop > PI/2) {
			periodAmpCrop = PI - periodAmpCrop;
		}
		periodCrop = 2 * PI / periodAmpCrop;
		setResult("Period (sec)_crop", i-1, periodCrop);
		setResult("Fitting (R^2)_crop", i-1, RCrop);
		close(); // For straighten line

		// Suggested oscillation mode
		setResult("Period ratio", i-1, period/periodCrop);
		decision = "undefined";
		if (period < maxPeriod && period > minPeriod) {
			if (mainDirection > minAngle || mainDirection < -minAngle) {
				decision = "travelling";
			} else if (period/periodCrop > minRatioP2P && period/periodCrop < maxRatioP2P) {
				decision = "pole-to-pole";
			} else if (period/periodCrop > minRatioPuls && period/periodCrop < maxRatioPuls) {
				decision = "pulsing";
			}
		}
		setResult("Suggested oscillation mode", i-1, decision);
	}
	
	if (answer[4] == 1) {
		close(); // For calibrated stack
	}
}

close(); // For original image
setBatchMode(false);

if (answer[3] == 1) {
	saveAs("Results", dir+title+".csv");
}
if (answer[2] == 0) {
	close("Results");
}