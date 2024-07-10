#Steps taken to merge SCRAM Movement Models as subdirectory into SCRAM2 repository

#Want to import SCRAM-Movement_Models as a subdirectory into SCRAM2 (subtree merge)
#if I were Andrew I could merge the two existing repositories together
#(assuming he creates or clones SCRAM-Movement-Models), but since I'm not Andrew, I need to clone SCRAM 2
# browseURL("https://happygitwithr.com/existing-github-last")
# browseURL("https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github")
# browseURL("https://www.theserverside.com/blog/Coffee-Talk-Java-News-Stories-and-Opinions/How-to-push-an-existing-project-to-GitHub")
# browseURL("https://gist.github.com/martinbuberl/b58fd967f271f32f51f50aee62e7332c")

## set directories
#dir.home <- file.path('P:') #Home directory (network)
#dir.root <- file.path(dir.home, 'SCRAM', 'BRI') #Root directory (network)

usr.name <- 'holly.goyert' #hard-coded (needs manual update)!
dir.home <- file.path('C:', 'Users', usr.name, 'OneDrive - Biodiversity Research Institute', 'Documents') #Home directory (local)
dir.root <- file.path(dir.home, 'ProjectsOffshoreWind', 'Atlantic', 'USFWS', 'SCRAM', 'BRI') #Root directory (local)

dir.create(file.path(dir.root, 'SCRAM2_shiny')) #subfolder for output files
dir.proj <- file.path(dir.root, 'SCRAM2_shiny') #Project directory

dir.work <- dir.proj
setwd(dir.work)

## Clone SCRAM 2 to local directory (need to start with existing repository) ####
#### need to start with existing repository, and 
#### don't want to overwrite server directory at P:\SCRAM\BRI\SCRAM2_shiny)

#browseURL("https://happygitwithr.com/new-github-first")
## 15 New project, GitHub first ####

## 15.2 New RStudio Project via git clone
usethis::create_from_github(
  "https://github.com/Biodiversity-Research-Institute/SCRAM2",
  destdir = dir.work
)
#Error: forking is disabled.
## Creates a new local directory in destdir, which is all of these things:
### a directory or folder on your computer
### a Git repository, linked to a remote GitHub repository
### an RStudio Project

## 15.2.2 RStudio IDE
## File > New Project > Version Control > Git
## Behind the scenes, RStudio has done this for you:
#git clone https://github.com/Biodiversity-Research-Institute/SCRAM2
## RStudio Project included (from clone)

# #in the terminal/shell (Tools>Terminal>New Terminal Alt+Shift+R):
# pwd
# cd SCRAM2 #Make this repo your working directory (cd is the command to change directory)
# ls #list its files

## Manually copy project files into the folder created by the clone.
## Perform a git add . and a git commit.
# git pull
# git add 'SCRAM Movement Models'/
# git commit -m "Add SCRAM 2.0 movement model scripts and data to repository for BOEM Report Addendum"
## Push your changes up to GitHub.
# git push

## 7/2/24 Opened: 
# C:\Users\holly.goyert\OneDrive - Biodiversity Research Institute\Documents\ProjectsOffshoreWind\Atlantic\USFWS\SCRAM\BRI\SCRAM2_shiny\SCRAM2\SCRAM2_shiny.Rproj
# P:\ #also updated (files copied over)
# browseURL("https://docs.github.com/en/repositories/working-with-files/managing-files/renaming-a-file")
## from within terminal:
# git mv 'SCRAM Movement Models/movement_model_report_addendum_feb2024_rekn.R' 'SCRAM Movement Models/movement_model_SCRAM2_1state_rekn.R'
# git mv 'SCRAM Movement Models/movement_model_report_addendum_feb2024_pipl.R' 'SCRAM Movement Models/movement_model_SCRAM2_1state_pipl.R'
# git mv 'SCRAM Movement Models/movement_model_report_addendum_feb2024_rost.R' 'SCRAM Movement Models/movement_model_SCRAM2_1state_rost.R'
# git mv 'SCRAM Movement Models/movement_model_report_addendum_feb2024_documentation.txt' 'SCRAM Movement Models/movement_model_report_fall2024_documentation.txt'

# git status
# git commit -m "SCRAM v2.1.4 rename files"
# git push origin main

# copied over from:
# BRI/Integrated Movement Models/ to:
# BRI/SCRAM2_shiny/SCRAM2/SCRAM Movement Models/
# 1. movement_model_report_fall2024_documentation.txt (revised from BRI/SCRAM2_shiny/SCRAM2/SCRAM Movement Models/movement_model_report_addendum_feb2024_documentation.txt)
# 2. serial_time.R
# 3. movement_model_SCRAM2_1state_rekn.R
# 4. movement_model_SCRAM2_1state_pipl.R
# 5. movement_model_SCRAM2_1state_rost.R

# pushed the following commit:
# SCRAM v2.1.4 calculates, fits, and predicts to 24-hour time step (rather than daily). 
# Replaces filtering functions with faster code to filter out false positive runs by site. 
# Fixed number of days in January for hourly data in serial_time.R

## Final model results for SCRAM v2.1.4 are contained within:
# C:\Users\holly.goyert\OneDrive - Biodiversity Research Institute\Documents\ProjectsOffshoreWind\Atlantic\USFWS\ and
# P:\
# ~\SCRAM\BRI\Integrated Movement Models\Motus Movement Models\REKN\REKN_2023_2021_1state_0D_var_hr24
# ~\SCRAM\BRI\Integrated Movement Models\Motus Movement Models\PIPL\PIPL_2023_2018_1state_0D_var_hr24
# ~\SCRAM\BRI\Integrated Movement Models\Motus Movement Models\ROST\ROST_2023_2018_1state_0D_var_hr24

# git mv 'SCRAM Movement Models/' 'Motus Movement Models/'
# git status
# git commit -m "SCRAM v2.1.4 rename subfolder"
# git push origin main

#Manually created 'GPS Movement Models' subfolder and copied over
# movement_model_SCRAM2_crawl_rekn_gps_out.R #(with slight edits for consistency with user.name from other scripts)
