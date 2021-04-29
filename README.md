# PFTC5-Intraspecific

## Downloading and 'cleaning' the data
This can be done by running the scripts/0_data_import.R file. This will then download the data as well as filter/clean the data so that we only have the data that we will be using for our analyses (i.e. 2020, sites (WAY, ACJ, TRE), control only sites, and our 6 species).

You can 'run'/call this script by calling source(here::here(path = "scripts/0_data_import.R")) at the beginning of your script. which will also add the two cleaned/filtered datasets to your environment as traits (for functional traits) as well as ensuring that you have the most up to date dataset stored on you machine.

Note if you do call source(here::here(path = "scripts/0_data_import.R")) at the beginning of your script it does mean you will always be using the most complete/updated dataset if we do end up changing which data we want to keep or exclude from analyses.


## Script naming structure
Each subtask related to/needing a coding workflow should be contained within its own script and should be named starting with the subtask number and a brief descriptive name. Scripts should be placed in the scripts/ folder. This means we can keep track of each task separately.

## Working on subtasks - using pull requests
Ideally each subtask should be on a new branch. This means that each subtask can be turned into a pull request (PR) allowing us to easily see the full commit history for that subtask and also allows subgroup members to request reviews/feedback from each other as well as have conversation threads. PRs can
initially be marked as drafts and once ready (i.e. completed) it can be marked as ready for review and then merged into the master branch.

branches should be named after the subtask code - same for the PR (although this can be a bit more comprehensive/descriptive).
