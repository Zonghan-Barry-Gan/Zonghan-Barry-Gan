macro "cell in gel"{

file='0';
//feed in root address    

root_address= "D:\\image_processing\\rename\\x10\\";
for (i=1;i<41;i++){
  //i=2
    address=root_address+i+"\\Pos0\\";
    address_1="D:\\image_processing\\newdo\\x10\\";
                print(address);
//open image files
    list1= getFileList(address);
//green    
    m=0;
    open(address+list1[m]);
    //run("Subtract Background...", "rolling=50 sliding");
    run("Enhance Contrast", "saturated pixels=5%");
    titleg="x10_"+i+"_0g.tiff";
    file_titleg=address_1+titleg;
    saveAs("Tiff",file_titleg);
    //close();
//red 
    m=1;
    open(address+list1[m]);
    //run("Subtract Background...", "rolling=50 sliding");
    run("Bandpass Filter...","filter large structures down to=40 pixels suppress stripes=vertical tolerance of direction=5% autoscale after filtering saturate image when autoscaling");
    titler="x10_"+i+"_1r.tiff";
    file_titler=address_1+titler;
    saveAs("Tiff",file_titler);
    //close();
//blue
    m=2;
    open(address+list1[m]);
    //run("Subtract Background...", "rolling=50 sliding");
    run("Bandpass Filter...","filter large structures down to=40 pixels filter small structures up to=3 pixels suppress stripes=vertical tolerance of direction=5% autoscale after filtering saturate image when autoscaling");
    titleb="x10_"+i+"_2b.tiff";
    file_titleb=address_1+titleb;
    saveAs("Tiff",file_titleb);
    //close();
//merge
run("Merge Channels...", "c1="  + titler + " c2="+ titleg + " c3=" + titleb + " keep");         
close(titler);
close(titleg);
close(titleb);
//run("RGB Color"); //Equivalent to Image > Type > RGB Color to fix the scale bar color
//x10 add a scale bar of 100um
run("Set Scale...", "distance=128.0625 known=100 pixel=1 unit=um global");
run("Scale Bar...", "width=100 height=8 font=28 color=White background=None location=[Lower Right] bold");       

//set scale and dpi to default 600-1.67-1.67
file_titlefinal=address_1+"final\\"+"x10_"+i+".tiff";
saveAs("Tiff",file_titlefinal);
//run("Scale to DPI", "create new window");
run("Close All");
    }
    }