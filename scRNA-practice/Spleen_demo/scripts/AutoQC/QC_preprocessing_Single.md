Notes:

The preprocessing steps include: 1)Soup removel using SoupX, 2)low quality cell removal using miQC, and 3) doublet removal using DoubletFinder.

The scripts are designed to run single samples. To submit the tasks, run the wrapper script by "./QC_preprocessing_wrapper.sh input/dir output/dir". The sample folders will be detected within the script and provided to the "QC_preprocessing_Single.sh" as sample names. The sample specific output directories will be created within the "Single" shell script, which functions as the job submission script.

The input directory should be the directory under which the samples are listed, and the alignment files generated (the "outs" directory) should be located under the sample directory.   
