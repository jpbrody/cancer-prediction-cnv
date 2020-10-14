Short version for experienced and impatient R users:
	1. Change the "project" varaible to your Google cloud project ID.
	2. Run AUC-ROC-plots.R   
        3. It uses your Google identity to download TCGA data from bigquery, so you need to authenticate yourself to Google.  A browser window will pop up asking for password/permission.  This query is about 1/4000th of the free monthly allocation provided by Google, so it won't cost you anything.
	4. The AUC-ROC-plots.R will generate ROC graphs and compute AUC values.

==================================================================
Long version:

The R script can be used to reproduce the results shown in the manuscript.

You need to first install R to run this script .  R is available from https://www.r-project.org/  This script was tested with version 3.5.0, but should work with any recent version.


The package H2O occasionally causes problems.  H2O requires java to run properly, and it doesn't work with java 9, the latest java release.  It needs java8.  See http://docs.h2o.ai/h2o/latest-stable/h2o-docs/faq/java.html

Once R and the packages are properly installed, you need to change a single line in the program.
You needs to be change 'project' to your Google cloud project ID.  See the code/comments below.

# This variable 'project' is the Google cloud project ID.  
project<- "test-project-cgc"
# You need to change it to *your* Google cloud project ID
# if you don't have one, go to http://bigquery.cloud.google.com/ and make one
# or see this page for help https://support.google.com/cloud/answer/6158840?hl=en
# The project ID is not secret.  You will need to authenticate yourself to google with your secret password.
# Google provides 1 TB of free queries per month.  This program only has one query that will consume less than 1 GB,
# so, it will use less than 1/1000th of your monthly quota.

# the first time you execute a Bigquery, a browser window will open and ask you to authenticate yourself with Google. 
# the authentication token will be cached, so you don't need to do it every time.


Installation time is minimal.  If you need to install R, it could take 5 minutes or so.
Run time is about 5-10 minutes on a 3 year old desktop.  There are progress indicators, so you should see graphs being generated, messages written to the console, etc.

