## Script reads the two NIST libraries and loads into sql
<?php
ini_set('memory_limit', '4096M');



$libdir = $argv[1];  // dir with the libs
$libraries = array('nist_msms','nist_msms2');



$link = mysqli_connect( 'localhost', '####', '####', '#####' );
$sql = "IF OBJECT_ID('nistlib', 'U') IS NOT NULL DROP TABLE nistlib"; 
mysqli_query($link, $sql);


$sql = "CREATE TABLE IF NOT EXISTS nistlib (".
		"library VARCHAR(50),".
		"LIBID TEXT,".
		"name VARCHAR(500),".
		"formula VARCHAR(100),".
		"MW DOUBLE(12,6),".
		"precursorMZ DOUBLE (12,6),".
		"casNO INT(10),".
		"ion VARCHAR(50),".
		"nPeaks INT(4),".
		"mzPeaks TEXT, ".
		"intPeaks TEXT )";
				
	#echo $sql . "\n";			
mysqli_query($link, $sql);


$i=0;
foreach($libraries as $libraryName){
	error_reporting(1);
	$sql="";
	
	echo("Reading nist_msms lib.....");
	echo($libdir . '\\' . $libraries[$i] . ".MSP\n\n\n");
	$data = file_get_contents($libdir . '\\' . $libraries[$i] . ".MSP");
    $arr = explode("Name: ",$data);
  
    echo("Patience...inserting library....\n");
	$j=0;
    foreach($arr as $chunk){

	echo($j."\r");
	$entry       = explode("\r\n",$chunk); #print_r($entry);
	if(count($entry)>1){
		     #print_r($entry);
	$name        = $entry[0];
	$LIBID       = intval(ltrim(explode(":",array_values(preg_grep('/^NISTNO/',      $entry))[0])[1]));
	
	if(strlen($LIBID) == 0){
		$LIBID = $entry[0];
	}
	#var_dump($LIBID);# . "\n";
	$formula     = ltrim(explode(":",array_values(preg_grep('/^Formula/',     $entry))[0])[1]);
	$MW          = ltrim(explode(":",array_values(preg_grep('/^MW/',          $entry))[0])[1]);
	$precursorMZ = ltrim(explode(":",array_values(preg_grep('/^PrecursorMZ/', $entry))[0])[1]);
	$casNO       = ltrim(explode(":",array_values(preg_grep('/^CASNO/',       $entry))[0])[1]);
	#$ion         =       explode("\$:03",array_values(preg_grep('/^Synon: \$:03/', $entry))[0])[1];
	$ion         = ltrim(explode(":",array_values(preg_grep('/^Precursor_type/',$entry))[0])[1]);
	$nPeaks      = ltrim(explode(":",array_values(preg_grep('/^Num peaks/',   $entry))[0])[1]);
	
	$IDX = array_keys(preg_grep('/^Num peaks/',  $entry))[0];
	$STARTIDX = $IDX +1;
	$peakArr = array_slice($entry,$STARTIDX,$nPeaks);
	
	$intArr = array();
	$mzArr = array();

	foreach($peakArr as $peak){
		$values = explode(" ",$peak);
		array_push($mzArr, $values[0]);
		array_push($intArr,$values[1]);
	}
	
	$mzPeaks = implode(" ",$mzArr);
	$intPeaks = implode(" ", $intArr);

	
	
	$sql = "INSERT into NISTLIB ( library,LIBID,   name,   formula,   MW,   precursorMZ,   casNO,   ion,   nPeaks, mzPeaks, intPeaks) ".
						"VALUES('$libraryName','$LIBID',\"$name\",'$formula','$MW','$precursorMZ','$casNO','$ion','$nPeaks','$mzPeaks','$intPeaks');\n";
    error_reporting(1);
	mysqli_query($link, $sql);
	$j++;
	}
   }

   echo("Now hold on there pal....just about done...\nNow Moving to next library. \n");
   

 
   $i++;
}

 mysqli_close($conn);
?>
