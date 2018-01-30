<?php
#### RUNs NIST's MSSearch program (as MSPEPSearch.exe):
### INSTALL MSPEPSEARCH FROM :  http://chemdata.nist.gov/dokuwiki/doku.php?id=peptidew:mspepsearch
###


#### Path to the BINARY
$BIN_PATH = "\\NIST17\\MSSEARCH\\mspepsearch.exe";

#### Path to the libraries......flat file in .msp (NIST) format
$LIB_PATH = "\\NIST17\\MSSEARCH";                                                                // Edited 4.19.2016  from nist to lib

#### THE WORKING DIRECTORY....stays the same throughout experiment/extraction
####  This is the location of the mgf files
$wd = $argv[1];

echo $DIR . "\n";

#### Get the mgf file list
$mgfFileList = glob( $wd . "\\*.mgf" );       								
print_r($mgfFileList);


##### Perform NIST search on each mgf!
#####   NIST WILL output searches in tsv format (tab delim)
#####   Please do not alter the search parameters and order of options.
#####   
foreach ( $mgfFileList as $MGF_FILE ) {
	
    # Remove nasty carriage return from file name.
	$fh = rtrim ( $MGF_FILE, ".mgf" );
	
	$OUTFILE = $fh . ".MGF.TSV\"";
	
	$LIB1 = "\"" . $LIB_PATH . "\\nist_msms\"";
	$LIB2 = "\"" . $LIB_PATH . "\\LipidBlast-pos\"";
	$LIB3 = "\"" . $LIB_PATH . "\\nist_msms2\"";
	
	$OPTIONS = " Gdi /OutSpecNum /OutPrecursorMZ /NumCompared /HiPri /MinMF 10 /PROGRESS /HITS 10 /Z 1.0 /M 0.8 /All ";
	$OPTIONS .= "/LIB $LIB1  /LIB $LIB2 /LIB $LIB3 ";
	$OPTIONS .= "/INP \"$MGF_FILE\" /OUTTAB \"$OUTFILE /OutNISTrn /OutPrecursorType /OutPrecursorMZ /OutSpecNum /OutCE /OutChemForm /OutCASrn";
	
	$CMD = $BIN_PATH . $OPTIONS;
	echo $CMD . "\n";
	
	system ( $CMD );
	
	
	
}



?>
