library(workflowr)
library(rmarkdown)
setwd(dir="/projects/compsci/vmp/USERS/heh/KOMP2022/")

# Configure Git (only need to do once per computer)
# wflow_git_config(user.name = "xhyuo", user.email = "hehao1986@gmail.com")
# wflow_use_github("KumarLabJax")

# Build the site
wflow_publish(c("analysis/*",
                "code/*",
                "data/*"
                ),
                "First build")

# publish data
wflow_publish("data/*.csv", "data",
              republish = FALSE)

wflow_status()

# git config --global user.name "xhyuo"
# git config --global user.email "hehao1986@gmail.com"
# git config -l
# git remote set-url origin https://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@github.com/KumarLabJax/QTL_analysis_for_Crichton.git

