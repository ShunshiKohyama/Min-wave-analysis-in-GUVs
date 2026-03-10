/**
 * Min_wave_analysis_in_GUVs (Fiji JavaScript Version)
 * * Description: 
 * Automated GUV detection and kymograph profile extraction for Min protein wave analysis.
 * This script uses a "Hybrid Architecture": JavaScript for robust GUI and dimension management, 
 * and an integrated IJM (ImageJ Macro) engine for high-speed, headless circle fitting.
 * Validated only for Fiji v1.54p or later
 * Input: 2D or 3D-stacked time-lapse image files with fluorescently labeled (inside or membrane) GUVs
 * Resulting data sets are saved as CSV files or a Zip file (for ROIs)
 * * Reference: Kohyama et al., Nature Communications, 2022.
 * * Author: Shunshi Kohyama (The University of Tokyo), 2026
 */

var IJ = Packages.ij.IJ;
var ImagePlus = Packages.ij.ImagePlus;
var WindowManager = Packages.ij.WindowManager;
var GenericDialog = Packages.ij.gui.GenericDialog;
var RoiManager = Packages.ij.plugin.frame.RoiManager;
var ResultsTable = Packages.ij.measure.ResultsTable;
var Duplicator = Packages.ij.plugin.Duplicator;
var OvalRoi = Packages.ij.gui.OvalRoi;
var WaitForUserDialog = Packages.ij.gui.WaitForUserDialog;
var DirectoryChooser = Packages.ij.io.DirectoryChooser;

// --- Initialize Metadata and Parameters ---
var imp = IJ.getImage();
var title = imp.getTitle();
var cal = imp.getCalibration();
var width = cal.pixelWidth;
var interval = cal.frameInterval;

var dims = imp.getDimensions(); // [width, height, nChannels, nSlices, nFrames]
var nC = dims[2]; var nZ = dims[3]; var nT = dims[4];

var methodsArray = Packages.ij.process.AutoThresholder.getMethods();

// --- Parameter Input ---
var gd = new GenericDialog("Min Wave Analysis Setup");
gd.addNumericField("Set scale (µm per pixel):", width, 4);
gd.addNumericField("Set interval (sec):", interval, 2);
gd.addNumericField("Set channel for analysis:", 1, 0);
gd.addNumericField("Set slice (z-position):", 1, 0);
gd.addNumericField("Set frame for detection:", 1, 0);
gd.addChoice("Threshold method:", methodsArray, "Huang");
gd.addNumericField("Min diameter (µm):", 15, 2);
gd.addNumericField("Max diameter (µm):", 50, 2);
gd.addNumericField("Min circularity:", 0.60, 2);
gd.addNumericField("Max circularity:", 1.00, 2);
gd.showDialog();

if (!gd.wasCanceled()) {
    var scale = gd.getNextNumber();
    interval = gd.getNextNumber();
    var channel = Math.min(gd.getNextNumber(), nC);
    var slice = Math.min(gd.getNextNumber(), nZ);
    var frame = Math.min(gd.getNextNumber(), nT);
    var type = gd.getNextChoice();
    var minDiam = gd.getNextNumber();
    var maxDiam = gd.getNextNumber();
    var minCirc = gd.getNextNumber();
    var maxCirc = gd.getNextNumber();
    
    cal.pixelWidth = scale; cal.pixelHeight = scale;
    cal.setUnit("µm");
    imp.setCalibration(cal);

    var rm = RoiManager.getInstance() || new RoiManager();
    rm.reset();
    
    var dupImp = new Duplicator().run(imp, channel, channel, slice, slice, frame, frame);
    var dupID = dupImp.getID();
    var satisfied = false; var skip = false; var range = 5;

    // --- Iterative GUV Detection ---
    while (!satisfied) {
        var minArea = Math.PI * Math.pow((minDiam / 2), 2);
        var maxArea = Math.PI * Math.pow((maxDiam / 2), 2);

        var maskIp = dupImp.getProcessor().duplicate();
        var maskImp = new ImagePlus("mask", maskIp);
        IJ.setAutoThreshold(maskImp, type + " dark");
        IJ.run(maskImp, "Convert to Mask", "");
        rm.reset();
        IJ.run(maskImp, "Analyze Particles...", "size=" + minArea + "-" + maxArea + " circularity=" + minCirc + "-" + maxCirc + " clear include add");
        maskImp.close();
        
        dupImp.show(); 
        IJ.selectWindow(dupID);
        rm.runCommand("Show All with labels");
        
        IJ.beep();
        var gd2 = new GenericDialog("Detection Review");
        gd2.addCheckbox("Proceed to circle fitting", true);
        gd2.addNumericField("Search range for fitting (µm):", 5, 2);
        gd2.addCheckbox("Adjust ROIs manually before fitting", false);
        gd2.addChoice("New Threshold method:", methodsArray, type);
        gd2.addNumericField("New Min diameter (µm):", minDiam);
        gd2.addNumericField("New Max diameter (µm):", maxDiam);
        gd2.addNumericField("New Min circularity:", minCirc);
        gd2.addNumericField("New Max circularity:", maxCirc);
        gd2.showDialog();

        if (gd2.wasCanceled()) break;
        satisfied = gd2.getNextBoolean();
        range = gd2.getNextNumber();
        skip = gd2.getNextBoolean();
        
        if (!satisfied) {
            type = gd2.getNextChoice();
            minDiam = gd2.getNextNumber();
            maxDiam = gd2.getNextNumber();
            minCirc = gd2.getNextNumber();
            maxCirc = gd2.getNextNumber();
            rm.runCommand("Remove Overlay");
            if (dupImp.getWindow()) dupImp.getWindow().setVisible(false);
        }
    }

    // Manual ROI adjustment if requested
    if (skip) {
        if (dupImp.getWindow()) dupImp.getWindow().setVisible(true);
        IJ.selectWindow(dupID);
        rm.runCommand("Show All with labels");
        new WaitForUserDialog("Manual ROI Adjustment", "Add or delete ROIs in the ROI Manager.\nClick OK when finished.").show();
    }

    // --- Native IJM Fitting Engine ---
    if (dupImp.getWindow()) dupImp.getWindow().setVisible(false);
    IJ.runMacro("setBatchMode(true);");
    IJ.showStatus("Executing high-speed circle fitting...");
    
    var pixels = Math.ceil(range / scale);
    var macroCode = 
        "setBatchMode(true);" + 
        "selectImage(" + dupID + ");" + 
        "roiManager('Show None');" + 
        "N = roiManager('count');" +
        "for (i=0; i<N; i++) {" +
        "  roiManager('select', 0);" +
        "  getSelectionBounds(x, y, w, h);" +
        "  d = minOf(w, h);" +
        "  pixels = " + pixels + ";" +
        "  if (pixels > d) pixels = d;" +
        "  roiManager('Delete');" +
        "  I = 0; X = x; Y = y; D = d;" +
        "  for (j=0; j<pixels; j++) {" +
        "    makeOval(x, y, d-j, d-j);" +
        "    run('Area to Line');" +
        "    for (k=0; k<=h-d+j; k++) {" +
        "      for (l=0; l<=w-d+j; l++) {" +
        "        Roi.move(x+l, y+k);" +
        "        getStatistics(area, mean);" +
        "        if (I < mean * area) { I = mean * area; X = x+l; Y = y+k; D = d-j; }" +
        "      }" +
        "    }" +
        "  }" +
        "  makeOval(X, Y, D, D); run('Area to Line'); roiManager('add');" +
        "}";
    
    IJ.runMacro(macroCode); 
    dupImp.changes = false; dupImp.close();

    // --- Metadata Export and Profile Extraction ---
    var dirChooser = new DirectoryChooser("Select Output Directory");
    var dir = dirChooser.getDirectory();
    if (dir) {
        var nameGD = new GenericDialog("Export Settings");
        nameGD.addStringField("Project Title Prefix:", title, 30);
        nameGD.showDialog();
        if (!nameGD.wasCanceled()) {
            title = nameGD.getNextString();

            // Export Metadata
            var rtMeta = new ResultsTable();
            var nROIs = rm.getCount();
            for (var i = 0; i < nROIs; i++) {
                rm.rename(i, title + "_" + (i + 1));
                rtMeta.incrementCounter();
                rtMeta.addValue("ID", rm.getRoi(i).getName() + ".csv");
                rtMeta.addValue("Diameter_um", rm.getRoi(i).getBounds().height * scale);
            }
            rtMeta.setValue("Scale_um_pixel", 0, scale);
            rtMeta.setValue("Interval_sec", 0, interval);
            rtMeta.setValue("Analyzed_Channel", 0, channel);
            rtMeta.setValue("Analyzed_Slice", 0, slice);
            rtMeta.save(dir + title + "_metadata.csv");

            var Straightener = Packages.ij.plugin.Straightener;
            var straightener = new Straightener();
            var nFramesActual = nT > 1 ? nT : nZ;
            var stack = imp.getStack();
            var dummyImp = new Packages.ij.ImagePlus("dummy", stack.getProcessor(1));

            for (var i = 0; i < nROIs; i++) {
                var rtProf = new ResultsTable();
                var roi = rm.getRoi(i);
                for (var j = 0; j < nFramesActual; j++) {
                    var stackIndex = (nT > 1) ? imp.getStackIndex(channel, slice, j + 1) : imp.getStackIndex(channel, j + 1, 1);
                    dummyImp.setProcessor(stack.getProcessor(stackIndex));
                    dummyImp.setRoi(roi);
                    var straightIp = straightener.straightenLine(dummyImp, 5); // 5-pixel width average
                    if (straightIp) {
                        var sW = straightIp.getWidth();
                        for (var x = 0; x < sW; x++) {
                            var sum = 0; for (var y = 0; y < 5; y++) sum += straightIp.getf(x, y);
                            rtProf.setValue(String(interval * j), x, sum / 5.0);
                        }
                    }
                }
                rtProf.save(dir + roi.getName() + ".csv");
            }
            rm.runCommand("Save", dir + title + "_ROI.zip");
        }
    }
    
    IJ.runMacro("setBatchMode(false);");
    IJ.selectWindow(imp.getID());
    rm.runCommand("Show All with labels");
    IJ.showStatus("Min Wave Analysis Complete.");
}
