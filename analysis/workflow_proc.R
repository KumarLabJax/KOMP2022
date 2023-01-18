library(workflowr)
library(rmarkdown)
setwd(dir="/projects/compsci/vmp/USERS/heh/KOMP2022/")

# Configure Git (only need to do once per computer)
# wflow_git_config(user.name = "xhyuo", user.email = "hehao1986@gmail.com")
# wflow_use_github("KumarLabJax")


# git config --global user.name "xhyuo"
# git config --global user.email "hehao1986@gmail.com"
# git config -l
# git remote set-url origin https://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@github.com/KumarLabJax/KOMP2022.git




# Build the site
wflow_git_commit(c("analysis/*",
                   "code/*",
                   "data/*"),
                "First build")

wflow_status()

# git config --global user.name "xhyuo"
# git config --global user.email "hehao1986@gmail.com"
# git config -l
# git remote set-url origin https://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@github.com/KumarLabJax/KOMP2022.git

