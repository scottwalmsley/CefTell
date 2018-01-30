<?php
###############################
### Script that reads all raw data (in .txt format)
##  from the working directory and then populates a new SQL table with this info.
## Used for grabbing raw spectra and plotting NIST search hits.



################################
###   raw MSMS filereader object
###    MAIN IS BELOW

class scanObj{
	private $index;
	private $scan;
	private $parentScan;
	private $mslevel;
	private $TIC;
	private $time;
	private $CeV;
	private $polarity;
	private $peakIntensity;
	private $isolationMZ;
	private $z;
	private $nPeaks;
	private $mzData;
	private $intData;

	public function getIndex() {
		return $this->index;
	}
	public function setIndex($index) {
		$this->index = $index;
		return $this;
	}
	public function getScan() {
		return $this->scan;
	}
	public function setScan($scan) {
		$this->scan = $scan;
		return $this;
	}
	public function getParentScan() {
		return $this->parentScan;
	}
	public function setParentScan($parentScan) {
		$this->parentScan = $parentScan;
		return $this;
	}
	public function getMslevel() {
		return $this->mslevel;
	}
	public function setMslevel($mslevel) {
		$this->mslevel = $mslevel;
		return $this;
	}
	public function getTIC() {
		return $this->TIC;
	}
	public function setTIC($TIC) {
		$this->TIC = $TIC;
		return $this;
	}
	public function getTime() {
		return $this->time;
	}
	public function setTime($time) {
		$this->time = $time;
		return $this;
	}
	public function getCeV() {
		return $this->CeV;
	}
	public function setCeV($CeV) {
		$this->CeV = $CeV;
		return $this;
	}
	public function getPolarity() {
		return $this->polarity;
	}
	public function setPolarity($polarity) {
		$this->polarity = $polarity;
		return $this;
	}
	public function getPeakIntensity() {
		return $this->peakIntensity;
	}
	public function setPeakIntensity($peakIntensity) {
		$this->peakIntensity = $peakIntensity;
		return $this;
	}
	public function getIsolationMZ() {
		return $this->isolationMZ;
	}
	public function setIsolationMZ($isolationMZ) {
		$this->isolationMZ = $isolationMZ;
		return $this;
	}
	public function getZ() {
		return $this->z;
	}
	public function setZ($z) {
		$this->z = $z;
		return $this;
	}
	public function getNPeaks() {
		return $this->nPeaks;
	}
	public function setNPeaks($nPeaks) {
		$this->nPeaks = $nPeaks;
		return $this;
	}
	public function getMzData() {
		return $this->mzData;
	}
	public function setMzData($mzData) {
		$this->mzData = $mzData;
		return $this;
	}
	public function getIntData() {
		return $this->intData;
	}
	public function setIntData($intData) {
		$this->intData = $intData;
		return $this;
	}

}

#################################################################
##### MAIN
#################################################################

### The working directory
$wd = $argv[1];

### The experimental tag:   typsically the name of the experiment....eg"30d_darpa"
$expTag = $argv[2];

### dbase
$dbase = $argv[3];

##################################################
### the main table for storing the RAW data from raw files
$link = mysqli_connect( 'localhost', 'root', 'camera9(', $dbase );
$query  = "CREATE TABLE IF NOT EXISTS raw_data_" ."$expTag (".

		"idx INT(8),".
		"scan INT(8),".
		"parentScan INT(8),".
		"mslevel INT(2),".
		"tic DOUBLE(13,4),".
		"scantime DOUBLE(12,4),".
		"cev INT(3),".
		"polarity VARCHAR(8),".
		"peakIntensity DOUBLE(12,2),".
		"isolationMZ DOUBLE(9,4),".
		"z INT(2),".
		"nPeaks INT(4),".
		"mzData TEXT,".
		"intData TEXT, ".
		"fh TEXT)";

mysqli_query($link, $query);


##### the raw ms txt files
$pwizfiles = glob($wd ."\\*.txt");

print_r($pwizfiles);#exit;

	  error_reporting(0);
foreach($pwizfiles as $fh){
	insertFileScanData($fh, $link,$expTag);
}

######   Function to extract file scan information from the raw .txt file and insert into tables.
function insertFileScanData($fh,$link,$expTag){
	$file = file_get_contents($fh);
	$arr = explode("      spectrum:",$file);

	foreach ( $arr as $chunk ) {

		$entry = array_map ( 'trim', explode ( "\n", $chunk ) );
		
		if (preg_match ( '/index:/', $entry [1] )) {

			$rawScanObj = new scanObj ();
			$rawScanObj->setIndex ( trim ( explode ( ":", $entry [1] )[1] ) );
			$rawScanObj->setScan ( trim ( explode ( "=", $entry [2] )[1] ) );
			$rawScanObj->setNPeaks ( trim ( explode ( ":", $entry [3] )[1] ) );
			$rawScanObj->setPolarity ( trim ( explode ( " ", $entry [4] )[1] ) );
			$rawScanObj->setTIC ( trim ( explode ( ",", $entry [7] )[1] ) );
			$rawScanObj->setTime ( trim ( explode ( ",", $entry [40] )[1] ) * 60.0 );
			if (preg_match ( '/cvParam: ms level, 1/', $entry [8] )) {
				#print_r($entry);
				$rawScanObj->setMslevel ( 1 );
				$rawScanObj->setMzData ( trim ( explode ( "]", $entry [47] )[1] ) );
				$rawScanObj->setIntData ( trim ( explode ( "]", $entry [50] )[1] ) );
			}
			if (preg_match ( '/cvParam: ms level, 2/', $entry [8] )) {
				$rawScanObj->setMslevel ( 2 );
				$rawScanObj->setCeV ( trim ( explode ( ",", $entry [57] )[1] ) );
				$rawScanObj->setPeakIntensity ( trim ( explode ( ",", $entry [54] )[1] ) );
				$rawScanObj->setParentScan ( trim ( explode ( "=", $entry [47] )[1] ) );
				$rawScanObj->setIsolationMZ ( trim ( explode ( ",", $entry [49] )[1] ) );
				$rawScanObj->setZ ( trim ( explode ( ",", $entry [53] )[1] ) );
				$rawScanObj->setMzData ( trim ( explode ( "]", $entry [60] )[1] ) );
				$rawScanObj->setIntData ( trim ( explode ( "]", $entry [63] )[1] ) );

				#echo trim( explode(",",$entry[71])[1]) . "\n";
				#print_r($entry);#exit;

			}
			$arr =  explode("\\",$fh);

		#print_r($arr);#exit;
		$filename = $arr[count($arr)-1];


		$query = "INSERT INTO raw_data_". "$expTag (" .
			"idx," .
			"scan," .
			"parentScan," .
			"mslevel," .
			"tic," .
			"scantime," .
			"cev," .
			"polarity," .
			"peakIntensity," .
			"isolationMZ," .
			"z," .
			"nPeaks," .
			"mzData," .
			"intData,".
			"fh) VALUES (" .
			"'".$rawScanObj->getIndex() ."'," .
			"'".$rawScanObj->getScan() ."'," .
			"'".$rawScanObj->getParentScan() ."'," .
			"'".$rawScanObj->getMslevel() ."'," .
			"'".$rawScanObj->getTIC() ."'," .
			"'".round($rawScanObj->getTime(),3) ."'," .
			"'".$rawScanObj->getCeV() ."'," .
			"'".$rawScanObj->getPolarity() ."'," .
			"'".$rawScanObj->getPeakIntensity() ."'," .
			"'".$rawScanObj->getIsolationMZ() ."'," .
			"'".$rawScanObj->getZ() ."'," .
			"'".$rawScanObj->getNPeaks() ."'," .
			"'".$rawScanObj->getMzData() ."'," .
			"'".$rawScanObj->getIntData() ."',".
			"'$filename')";
			#echo $query . "\n";
			mysqli_query($link, $query);
		}

	}

}


?>
