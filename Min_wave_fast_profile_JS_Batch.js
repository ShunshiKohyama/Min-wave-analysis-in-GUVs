/**
 * Min_wave_analysis_in_GUVs (Fiji JavaScript Version)
 * * Description: 
 * Batch processing of the automated GUV detection and kymograph profile extraction for Min protein wave analysis.
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
var GenericDialog = Packages.ij.gui.GenericDialog;
var RoiManager = Packages.ij.plugin.frame.RoiManager;
var ResultsTable = Packages.ij.measure.ResultsTable;
var Duplicator = Packages.ij.plugin.Duplicator;
var OvalRoi = Packages.ij.gui.OvalRoi;
var DirectoryChooser = Packages.ij.io.DirectoryChooser;
var Font = java.awt.Font;
var File = java.io.File;

// --- Choose directry ---
var dirChooser = new DirectoryChooser("Select Folder with TIF files");
var dir = dirChooser.getDirectory();
if (!dir) {
    IJ.log("Batch processing canceled.");
} else {
    var folder = new File(dir);
    var listOfFiles = folder.listFiles();
    var tifFiles = [];
    
    for (var i = 0; i < listOfFiles.length; i++) {
        if (listOfFiles[i].isFile() && listOfFiles[i].getName().toLowerCase().endsWith(".tif")) {
            tifFiles.push(listOfFiles[i]);
        }
    }

    if (tifFiles.length === 0) {
        IJ.showMessage("Error", "No TIF files found in the selected directory.");
    } else {
        var firstImp = IJ.openImage(tifFiles[0].getAbsolutePath());
        var width = 1.0; var interval = 1.0; var unit = "pixel";
        var nC = 1, nZ = 1, nT = 1;
        
        if (firstImp) {
            var cal = firstImp.getCalibration();
            width = cal.pixelWidth;
            interval = cal.frameInterval;
            unit = cal.getUnit();
            var dims = firstImp.getDimensions();
            nC = dims[2]; nZ = dims[3]; nT = dims[4];
            firstImp.close();
        }

        var methodsArray = Packages.ij.process.AutoThresholder.getMethods();

        // --- Initialize Metadata and Parameters ---
        var gd = new GenericDialog("Batch Min Wave Analysis Setup");

        if (unit != "µm") gd.addMessage("                           Please check!! (Unit is not µm)", new Font("SansSerif", Font.BOLD, 12));
        gd.addNumericField("Set scale (µm per pixel):", width, 4);

        if (interval == 0 && nT > 1) gd.addMessage("                           Please check!! (Frame Interval is 0)", new Font("SansSerif", Font.BOLD, 12));
        gd.addNumericField("Set interval (sec):", interval, 2);

        gd.addNumericField("Set channel for analysis:", 1, 0);
        gd.addNumericField("Set slice (z-position):", 1, 0);
        gd.addNumericField("Set frame for detection:", 1, 0);
        gd.addChoice("Threshold method:", methodsArray, "Huang");
        gd.addNumericField("Min diameter (µm):", 15, 2);
        gd.addNumericField("Max diameter (µm):", 50, 2);
        gd.addNumericField("Min circularity:", 0.60, 2);
        gd.addNumericField("Max circularity:", 1.00, 2);
        gd.addNumericField("Search range for fitting (µm):", 5, 2); 
        gd.showDialog();

        if (!gd.wasCanceled()) {
            var scale = gd.getNextNumber();
            var finalInterval = Math.round(gd.getNextNumber() * 100) / 100;
            var channel = gd.getNextNumber();
            var slice = gd.getNextNumber();
            var frame = gd.getNextNumber();
            var type = gd.getNextChoice();
            var minDiam = gd.getNextNumber();
            var maxDiam = gd.getNextNumber();
            var minCirc = gd.getNextNumber();
            var maxCirc = gd.getNextNumber();
            var range = gd.getNextNumber();
            var pixels = Math.ceil(range / scale);

            var rm = RoiManager.getInstance() || new RoiManager();

            // --- Batch Processing Loop ---
            for (var f = 0; f < tifFiles.length; f++) {
                var file = tifFiles[f];
                var fileName = file.getName();
                var baseName = fileName.substring(0, fileName.lastIndexOf(".")); 
                
                IJ.showStatus("Processing " + (f + 1) + "/" + tifFiles.length + ": " + fileName);
                
                var imp = IJ.openImage(file.getAbsolutePath());
                if (!imp) continue;

                var impCal = imp.getCalibration();
                impCal.pixelWidth = scale; impCal.pixelHeight = scale;
                impCal.setUnit("µm");
                imp.setCalibration(impCal);
                
                rm.reset();

                var curChannel = Math.min(channel, imp.getNChannels());
                var curSlice = Math.min(slice, imp.getNSlices());
                var curFrame = Math.min(frame, imp.getNFrames());

                var dupImp = new Duplicator().run(imp, curChannel, curChannel, curSlice, curSlice, curFrame, curFrame);
                
                dupImp.show(); 
                if (dupImp.getWindow()) dupImp.getWindow().setVisible(false);
                var dupID = dupImp.getID();

                // --- GUV Detection ---
                var minArea = Math.PI * Math.pow((minDiam / 2), 2);
                var maxArea = Math.PI * Math.pow((maxDiam / 2), 2);
                var maskIp = dupImp.getProcessor().duplicate();
                var maskImp = new ImagePlus("mask", maskIp);
                maskImp.setCalibration(dupImp.getCalibration());
                
                Packages.ij.Prefs.blackBackground = true; 
                IJ.setAutoThreshold(maskImp, type + " dark");
                IJ.run(maskImp, "Convert to Mask", "background=Dark"); 
                
                var analyzeCmd = "size=" + IJ.d2s(minArea,2) + "-" + IJ.d2s(maxArea,2) + " circularity=" + IJ.d2s(minCirc,2) + "-" + IJ.d2s(maxCirc,2) + " clear include add";
                IJ.run(maskImp, "Analyze Particles...", analyzeCmd);
                maskImp.close();

                // --- Native IJM Fitting ---
                var macroCode = 
                    "setBatchMode(true);" + 
                    "selectImage(" + dupID + ");" + 
                    "roiManager('Show None');" + 
                    "N = roiManager('count');" +
                    "for (i=0; i<N; i++) {" +
                    "  roiManager('select', 0);" +
                    "  getSelectionBounds(x, y, w, h);" +
                    "  d = minOf(w, h);" +
                    "  px = " + pixels + ";" +
                    "  if (px > d) px = d;" +
                    "  roiManager('Delete');" +
                    "  I = 0; X = x; Y = y; D = d;" +
                    "  for (j=0; j<px; j++) {" +
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
                var nROIs = rm.getCount();
                if (nROIs > 0) {
                    var rtMeta = new ResultsTable();
                    for (var i = 0; i < nROIs; i++) {
                        var padNum = (i + 1 < 10) ? "0" + (i + 1) : "" + (i + 1);
                        var roiName = baseName + "_" + padNum;
                        
                        rm.rename(i, roiName);
                        rtMeta.incrementCounter();
                        rtMeta.addValue("ID", roiName + ".csv");
                        rtMeta.addValue("Diameter_um", rm.getRoi(i).getBounds().height * scale);
                    }
                    rtMeta.setValue("Scale_um_pixel", 0, scale);
                    rtMeta.setValue("Interval_sec", 0, finalInterval);
                    rtMeta.setValue("Analyzed_Channel", 0, curChannel);
                    rtMeta.setValue("Analyzed_Slice", 0, curSlice);
                    rtMeta.save(dir + baseName + "_metadata.csv");

                    var Straightener = Packages.ij.plugin.Straightener;
                    var straightener = new Straightener();
                    var nFramesActual = imp.getNFrames() > 1 ? imp.getNFrames() : imp.getNSlices();
                    var stack = imp.getStack();
                    var dummyImp = new Packages.ij.ImagePlus("dummy", stack.getProcessor(1));

                    for (var i = 0; i < nROIs; i++) {
                        var rtProf = new ResultsTable();
                        var roi = rm.getRoi(i);
                        for (var j = 0; j < nFramesActual; j++) {
                            var stackIndex = (imp.getNFrames() > 1) ? imp.getStackIndex(curChannel, curSlice, j + 1) : imp.getStackIndex(curChannel, j + 1, 1);
                            dummyImp.setProcessor(stack.getProcessor(stackIndex));
                            dummyImp.setRoi(roi);
                            var straightIp = straightener.straightenLine(dummyImp, 5); 
                            
                            if (straightIp) {
                                var sW = straightIp.getWidth();
                                var currentTime = Math.round((finalInterval * j) * 100) / 100;
                                for (var x = 0; x < sW; x++) {
                                    var sum = 0; for (var y = 0; y < 5; y++) sum += straightIp.getf(x, y);
                                    rtProf.setValue(String(currentTime), x, sum / 5.0);
                                }
                            }
                        }
                        rtProf.save(dir + roi.getName() + ".csv");
                    }
                    rm.runCommand("Save", dir + baseName + "_ROI.zip");
                }
                
                imp.changes = false;
                imp.close();
            }

            rm.reset();
            IJ.showStatus("Batch Min Wave Analysis Complete.");
            IJ.showMessage("Finished", "Batch processing is complete!\nResults are saved in:\n" + dir);
        }
    }
}
