//This macro is used for quantifying intensity within the nuclear masks

image = getTitle();

roiManager("reset");

run("Clear Results");

run("Split Channels");

selectWindow("C1-"+image);

waitForUser("Please adjust brightness and contrast");

run("Median...", "radius=10");

setAutoThreshold("Mean dark");

run("Convert to Mask");

run("Watershed");

run("Analyze Particles...", "size=50-Infinity display exclude clear add");

selectWindow("C2-"+image);

run("Set Measurements...", "area integrated redirect=None decimal=3");

run("Clear Results");

roiManager("deselect");

roiManager("Measure");

waitForUser("Please copy data from the result table for nuclear intensity");

run("Close All");

print("finish analyzing "+image);