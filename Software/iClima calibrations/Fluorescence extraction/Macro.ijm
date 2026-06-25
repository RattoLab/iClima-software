close("*");
directory = getDir("Choose a Directory");
filelist = getFileList(directory);
//for (i = 0; i < lengthOf(filelist); i++) {
for (i = 0; i < lengthOf(filelist); i++) {
    if (endsWith(filelist[i], ".lif")) { 
        open(directory + File.separator + filelist[i]);
        out_dir = directory + File.separator + File.getNameWithoutExtension(filelist[i]) + File.separator ;
        File.makeDirectory(out_dir);
        n=nImages;
        for (i_image = 1;  i_image<=n; i_image++){
        	selectImage(i_image);
        	//nn = String.format("%03d",i_image);
        	title = getTitle() + "#" + i_image +"#";
        	print(title);
        	saveAs("tiff", out_dir+title);
        }
        close("*");
    } 
}
print("Done");