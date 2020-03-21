Run in the terminal:

module load matlab/R2015b

java -Xmx8000M -cp /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/BIODICA_GUI.jar BIODICAPipeLine -config /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/config -datatable mgen_caste_eset.txt  -doicamatlab 5,10,20,25,30,35,40,45

java -Xmx8000M -cp /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/BIODICA_GUI.jar BIODICAPipeLine -config /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/config -datatable mgen_caste_eset.txt  -donumbercomponents stability

java -Xmx8000M -cp /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/BIODICA_GUI.jar BIODICAPipeLine -config /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/config -datatable mgen_caste_abdO_eset.txt  -doicamatlab 5,10,15,20

java -Xmx8000M -cp /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/BIODICA_GUI.jar BIODICAPipeLine -config /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/config -datatable mgen_caste_abdO_eset.txt  -donumbercomponents stability

java -Xmx8000M -cp /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/BIODICA_GUI.jar BIODICAPipeLine -config /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/config -datatable mgen_dev_eset.txt  -doicamatlab 5,10,20,25,30,35,40,45

java -Xmx8000M -cp /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/BIODICA_GUI.jar BIODICAPipeLine -config /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/config -datatable mgen_dev_eset.txt  -donumbercomponents stability

java -Xmx8000M -cp /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/BIODICA_GUI.jar BIODICAPipeLine -config /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/config -datatable mgen_dev_eset.txt -outputfolder mgen_dev_eset_47IC -doicamatlab 47

java -Xmx8000M -cp /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/BIODICA_GUI.jar BIODICAPipeLine -config /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/config -datatable mgen_caste_eset.txt -outputfolder mgen_caste_eset_12IC -doicamatlab 12

java -Xmx8000M -cp /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/BIODICA_GUI.jar BIODICAPipeLine -config /uufs/chpc.utah.edu/sys/installdir/biodica/BIODICA/config -dobbhgraph bbh/
