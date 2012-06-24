afgwardiary
===========

Purpose:    Replication code for the paper "Point process modelling of the Afghan War Diary".

            The raw AWD is not available for download via this repository and only the minimum required 
            for replication (time and geo-tags only) are supplied. Required data resides in AfghanDataAllDay.mat
            (MATLAB), DatasetS1 (csv) and ./TextS5_Corroboration/origData (R software).
            
Usage:      Each folder contains MATLAB or R files which are self-contained scripts. When more than one file is present or the file to run is not obvious, dedicated readme files are present with further instructions. The three most important/relevant folders in this work are Main_VBAnalysis, Main_ResultsSection and Main_PredictionSection with code for the inference, results plotting and prediction respectively. TextS6_BasisSelection contains code for selecting the basis used after the PACF/PCCF analysis using code in TextS6_SpatialAnalysis. The names of the other folders are indicative of which section in the paper the code they contain is relevant to. 

Project page: http://andrewzm.github.com/afgwardiary/

Contact:    Andrew Zammit Mangion

Email:      azammit2 (at) inf.ed.ac.uk

Date:       2012-06-13


Copyright (c) 2012, under the Simplified BSD License. 

For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php
All rights reserved.