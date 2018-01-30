#Script loads NIST search results from a MSMS search to a SQL database

<?php
ini_set('memory_limit', '4096M');
###########################################
### point to your includes....theses are necesaary files.
include "sampleMetaObj.php";
include "experimentDbMetaSQL.php";

$dbase = $argv[6];            // Edited 4.22.2016 to add dbase to be passed to this program

#### hool up to the DB
$link = mysqli_connect( 'localhost', '####', '#######',$dbase);
#echo($link);
#exit;
##################################################
## Loads MSMS data from NIST MS search results (tsv files) into SQL sdb.
##################################################
#### the experiment tag   Use the same tag for all expts.  eg "30d_darpa_aqmeoh_pos"
$expTag = $argv[1];

### Extraction type shall be "AQ_MEOH","AQ_MTBE","LIPID_MTBE"....stays constant.
$extractionType = $argv[2];

### The polarity of the raw data files being parsed.   shall always be "pos" or "neg".
###   When analysing a batch of files, extraction type and polarities are kept separate, and separate table created.
###     eg...."AQ_MEOH and pos" files are analyzed separately from "AQ_MEOH" and"neg" files, and separate tables are created for theses.
$polarity = $argv[3];

### the working directory containing the nist results and raw files...etc....
$wd = $argv[4];

### true or false, create a table or not.   likely to be deprecated.
$create = $argv[5];

if($create == true){
$query = "CREATE TABLE IF NOT EXISTS ". $expTag ." (".
		"fh VARCHAR(300),".
		"extraction VARCHAR(20),".
		"scan INT(8),". //
		"polarity TINYTEXT,".
		"rt DOUBLE(8,3),". //
		"formula TINYTEXT,".
		"ion TINYTEXT,".
		"libz INT(2),".
		"libID TEXT,".
		"libMZ double(8,3),".
		"z INT(2),". //
		"pAb DOUBLE(14,2),". //
		"mz DOUBLE(8,3),".
		"cev INT(3),".
		"idRank INT(10),".
		"lib VARCHAR(30),".
		"score INT(3),".
		"dotP INT(4),".
		"prob DOUBLE(4,1),".
		"name VARCHAR(400))";

		mysqli_query($link, $query);
}
echo($query . "\n");#exit;
	$tsvfiles = glob ( $wd . "\\*.TSV" );
    print_r($tsvfiles);
	error_reporting(0);
	foreach ( $tsvfiles as $fh ) {

		#############################################################
		#open the pwiz text data file and get the information


		$pwfileName = str_replace(".MGF.TSV", ".txt", $fh);

		echo $pwfileName . "\n";
		$mzFile = file ( $pwfileName );
		echo $fh . "\n";

		list ( $end ) = array_keys ( array_slice ( $mzFile, - 1, 1, true ) );

		// ## Iterate the scanIDs and get relevent lines from the data file

		$ms2scanID = array ();
		$pAB = array ();
		$Cev = array ();

		$i = 0;

		while ( $i < count ( $mzFile ) ) {
			if (preg_match ( '/^id: scanId=/', ltrim ( $mzFile [$i] ) )) {
				#echo($mzFile[$i]);
				$ms2scanNum = trim ( explode ( "=", $mzFile [$i] )[1] );
                #cho($ms2scanNum);
				while ( ! preg_match ( '/^spectrum:/', ltrim ( $mzFile [$i] ) ) && ($i !== $end) ) {

					if (preg_match ( '/^cvParam: base peak intensity/', ltrim ( $mzFile [$i] ) )) {
						$Ab = trim ( explode ( ",", $mzFile [$i] )[1] );
					}

					if (preg_match ( '/^cvParam: ms level/', ltrim ( $mzFile [$i] ) )) {
						$level = trim ( explode ( ",", $mzFile [$i] )[1] );
						if ($level == 2) {
							while ( ! preg_match ( '/^binaryDataArray:/', ltrim ( $mzFile [$i] ) ) ) {

								if (preg_match ( '/^cvParam: collision energy/', ltrim ( $mzFile [$i] ) )) {

									$energy = trim ( explode ( ",", $mzFile [$i] )[1] );
									array_push ( $Cev, $energy );
									array_push ( $pAB, $Ab );
									array_push ( $ms2scanID, $ms2scanNum );
								}
								$i ++;
							}
						}
					}
					$i ++;
				}
			}
			$i ++;
		}


		########################################################################################
		##  Now get the search results from the nist file and link the cev and ab to the results
		########################################################################################
		$j = 0;
		$data = file($fh);
		$end = count($data)-6;
		$data = array_slice($data, 4, $end);
		echo "\n";

		foreach($data as $line){
		$la = explode("\t",$line);
		$la2 = explode(" ", $la[1]);
	    #print_r($la2);

		
		$scan = explode("=",$la2[0])[1];
		$rt = explode(":",$la2[1])[1];
		$z = $la[12];
        #echo("$rt $scan $z\n"); 
		$searchkey = array_search($scan,$ms2scanID);
		$energy = $Cev[$searchkey];
		$Ab = $pAB[$searchkey];


		$mz = $la[2];

		if(isset($la[3])){

			#print_r($la);
		    $idrank = $la[3];
			$score = $la[8];
			$dotP = $la[9];
			$prob = $la[10];
			$name = $la[11];
			$lib = $la[4];
			$libID = $la[20];
			$libz = $la[12];
			$libMZ = $la[7];

			$formula = $la[13];
			$ionType =  $la[14];

		}else{
		$idrank = $score = $prob = $name = $lib = $libz = $formula = $ionType = $libID = $libMZ = $dotP = NULL;
		}

		$tmp = explode("\\",$fh);
		$fileName = $tmp[count($tmp)-1];

		$query = "INSERT INTO $expTag (fh,extraction, polarity, scan,rt,formula,ion,libz,z,mz,pAb,cev,idRank,lib,libID,libMZ,score,dotP,prob,name) ".
		"VALUES ('$fileName','$extractionType','$polarity','$scan', '$rt','$formula','$ionType','$libz','$z','$mz','$Ab','$energy','$idrank',\"$lib\",'$libID','$libMZ','$score','$dotP','$prob',\"$name\")";
       # $query = "INSERT INTO $expTag (fh,extraction, polarity,scan,rt) VALUES ('$fileName','$extractionType','$polarity',$scan,$rt)";
		#echo $query . "\n\n";

		echo $j . "\r";
		$j++;
		mysqli_query($link, $query);
		#exit;

	}

	}
	error_reporting(1);


	?>



