dir1 = "/Users/oaksgil/Desktop/20250715_GO108_projections/";
list = getFileList(dir1);

// Output file in same directory
ofile = File.open(dir1 + "output_results_go108.txt");
print(ofile, "stack_name\tnum_roi\tdapi_int\tedu_int\tpla_ct\tab_int\n");

for (i = 0; i < list.length; i++) {
    if (!endsWith(list[i], ".dv")) continue;  // Skip non-DV files

    filename = dir1 + list[i];
    print("Processing file: " + filename);

    // Open all series in the file (i.e., all stacks)
    s = "open=[" + filename + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT open_all_series";
    run("Bio-Formats Importer", s);

    // Get titles of open windows (one per series)
    titles = getList("image.titles");
    for (k = 0; k < titles.length; k++) {
        selectImage(titles[k]);
        stackName = getTitle();
        print("Processing stack: " + stackName);

        // --- NEW: Duplicate channel 1 only for ROI creation ---
        Stack.setChannel(1);
        run("Duplicate...", "title=dup_" + stackName + " channels=1");
        dupTitle = getTitle(); 
        selectImage(dupTitle);

        run("8-bit");
        run("Gaussian Blur...", "sigma=1");
        run("Enhance Contrast...", "saturated=0.5");
        setAutoThreshold("Default dark");
        setOption("BlackBackground", true);
        run("Convert to Mask", "method=Shanbhag background=Dark calculate only black");
        run("Watershed", "slice");
        

        run("Analyze Particles...", "size=2500-10000 pixel circularity=0.8-1.00 display exclude clear add");
        
        selectImage(dupTitle);
       	close();
        

        num_roi = roiManager("count");

        for (j = 0; j < num_roi; j++) {
            roiManager("Select", j);
			
			selectImage(stackName);
            run("Clear Results");
            
            Stack.setChannel(1); // DAPI
            run("Measure");
            dapi_int = getResult("Mean", nResults - 1);
			
            Stack.setChannel(2); // EdU
            run("Measure");
            edu_int = getResult("Mean", nResults - 1);
			
            Stack.setChannel(3);  // PLA foci channel
            //run("Sharpen", "slice");
            run("Enhance Contrast...", "saturated=0.1");
            run("Find Maxima...", "prominence=800 strict output=[Count]");
            n = getResult("Count", nResults - 1);
			
            Stack.setChannel(4); // Antibody
            run("Measure");
            ab_int = getResult("Mean", nResults - 1);

            print(ofile, stackName + "\t" + j + "\t" + dapi_int + "\t" + edu_int + "\t" + n + "\t" + ab_int + "\n");
        }

        roiManager("reset");
        run("Clear Results");
        close(stackName);           // Close the original stack
    }
}

File.close(ofile);
print("Done.");