# KOMP2022

A [workflowr][] project.

[workflowr]: https://github.com/jdblischak/workflowr


The project folder is in the sumner /projects/compsci/vmp/USERS/heh/KOMP2022/.

raw data from LIMCS 2022-05-13.zip. I split the raw data into two zip files (2022-05-13-0.z01 and 2022-05-13-0.zip) less than 100mb using command:

zip 2022-05-13.zip --out 2022-05-13-0 -s 90m

use this command to unzip:

zip -F 2022-05-13-0.zip --out 2022-05-13.zip; unzip 2022-05-13.zip

Note: use the date name in the following Rscript for differentiating batches. for example here: batch = "2022-05-13"

If you run all the code in your local, please look over the parameter "projectDir and change it, in analysis/preprocessing_step1.R, analysis/preprocessing_step2.R and code/RCode_Functions_v3.R.

we are going to preprocess the raw data by running preprocessing_step1.R first and then running preprocessing_step2.R

cd /projects/compsci/vmp/USERS/heh/KOMP2022/analysis; nohup Rscript preprocessing_step1.R > preprocessing_step1.Rout 2>&1 &

cd /projects/compsci/vmp/USERS/heh/KOMP2022/analysis; nohup Rscript preprocessing_step2.R > preprocessing_step2.Rout 2>&1 &

Then open R file "analysis/run_AT_AllModels_script.R" and run it. it will generate the
"analysis/run_AT_AllModels_script.sh" file.

Submit the job by the following command.

"cd /projects/compsci/vmp/USERS/heh/KOMP2022/analysis; chmod a+x run_AT_AllModels_script.sh; ./run_AT_AllModels_script.sh"

KOMP_AT.Rmd will generate all the figures for the association results.
