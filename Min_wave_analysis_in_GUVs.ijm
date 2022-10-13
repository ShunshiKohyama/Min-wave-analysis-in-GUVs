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
Stack.getDimensions(widthPixel, heightPixel, channels, slices, frames);
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
close(); // Close the Mask
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
close(); // Close the duplicated image for detection

// Manually add or delete ROIs (here newly added ROIs must be circular lines but not circles by running the "Area to Line" command)
setBatchMode("show");
Overlay.hide;
roiManager("show all with labels");
beep();
waitForUser("Add or delete ROIs");

// Choose the actions and set searching range for floating calibration
Dialog.create("Analyze setting");
labels = newArray("Floating calibration", "Save calibrated stack", "Save kymographs", "Save ROIs");
defaults = newArray(false, false, true, false);
answer = newArray(4);
Dialog.addCheckboxGroup(2, 2, labels, defaults);
Dialog.addNumber("Searching range (floating calibration) (µm)", 2);
Dialog.addNumber("Channel for calibration", 1);
Dialog.show();
for (i=0; i<4; i++) {
	answer[i] = Dialog.getCheckbox();
}
range = Dialog.getNumber();
calibChannel = Dialog.getNumber();
n = nSlices;
N = roiManager("count");

if (answer[1] == 1 || answer[2] == 1 || answer[3] == 1) {
	dir = getDirectory("Choose a Directory");
}
	
if (answer[0] == 1) { // Make arrays for floating calibration
	X = newArray(n);
	Y = newArray(n);
	name = getTitle();
	Overlay.hide;
}

setBatchMode("hide");
run("Duplicate...", " ");
run("Clear Results");
if (answer[3] == 1) {
	roiManager("save", dir+title+"_ROI.zip") // Save the set of ROIs
}

for (i=1; i<=N; i++) {
	roiManager("select", i-1);

	// Floating calibration
	if (answer[0] == 1) {
		Stack.setChannel(calibChannel);
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
		
		if (answer[1] == 1) {
			saveAs("Tiff", dir+title+"_"+i+".tif");
		}
		if (answer[2] == 1) {
			makeOval(pixels, pixels, w, h); // In the case of making kymograph, make a circle along with the peripheral and then make a circular line
			run("Area to Line");
		}
	}

	// Make kymographs
	if (answer[2] == 1) {
		run("Straighten...", "title=straight line=5 process"); // Convert the circular line into the straight line with 5 pixel thickness along with the line
		w = getValue("Width"); // Get the length of the straightened line in pixel
		run("Size...", "width=w height=1 depth=n average interpolation=Bilinear"); // Convert the image size into 1 pixel height (with the same width)
		run("Make Montage...", "columns=1 rows=n scale=1"); // Make a kymograph by periperal intensity (width) vs time-points (height)
		if (answer[2] == 1) {
			saveAs("Tiff", dir+title+"_"+i+"_kymo.tif");
			if (answer[0] == 0) {
				close(); // For kymograph
				close(); // For straight line
			}
		}
	}
		
	if (answer[0] == 1) {
		close(); // For calibrated stack
	}
}

close(); // For original image
setBatchMode(false);