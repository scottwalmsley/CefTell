# Reads the lipid blast database and loads into SQL
<?php


$libdir = $argv[1];   //directory path to lipid-blast libs
$libraries = array("LipidBlast-pos");

#$link = mysqli_connect( 'RDLabData', 'root', 'J0ecamel', 'darpa_18' );
$link = mysqli_connect( 'localhost', 'root', 'camera9(', 'baldb' );


$i=0;
foreach($libraries as $libraryName){
	error_reporting(0);
	$sql="";
	
	echo("Reading nist_msms lib.....\n\n");
	echo($libdir . '\\' . $libraries[$i] . ".MSP\n\n\n");
	$data = file_get_contents($libdir . '\\' . $libraries[$i] . ".MSP");
    $arr = explode("Name: ",$data);
    
    echo("Patience...inserting library....\n");
	$j=0;
foreach($arr as $chunk){
	
	echo($j."\r");
	$entry       =  explode("\r\n",$chunk);
	#print_r($entry);
	$name        = ltrim(explode(";",$entry[0])[2]);
	$LM_LIBID    = $entry[0];

	$commentline     = array_values(preg_grep('/^Comment/',     $entry));
	$commentlineArr  = explode(";",$commentline[0]);
	#print_r($commentlineArr);
	$formula         = ltrim($commentlineArr[4]);
	$ion             = ltrim($commentlineArr[2]);
    $subArr          = explode(" ",$commentlineArr[0]);
	$MW = doubleval(explode("=", $subArr[2])[1]);
	$precursorMZ = doubleval(explode("=", $subArr[1])[1]);

	$casNO       = intval(ltrim(explode(":",array_values(preg_grep('/^CASNO/',       $entry))[0])[1]));
	$nPeaks      = intval(ltrim(explode(":",array_values(preg_grep('/^Num peaks/',   $entry))[0])[1]));

	#echo ("$name | $LM_LIBID || $formula || $ion || $MW || $precursorMZ || $casNO || $nPeaks\n");


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

	$sql = "INSERT into nistlib ( library, LIBID,   name , formula,MW,precursorMZ,casNO,ion,nPeaks,mzPeaks,intPeaks) VALUES ('$libraryName','$LM_LIBID','$name', '$formula','$MW','$precursorMZ','$casNO','$ion','$nPeaks','$mzPeaks','$intPeaks')";
	#$sql = "INSERT into nistlib ( library, LIBID,   name,   formula,   MW,   precursorMZ,   casNO,   ion,   nPeaks, mzPeaks, intPeaks) ".
	#					"VALUES('$libraryName','$LM_LIBID','$name','$formula','$MW','$precursorMZ','$casNO','$ion','$nPeaks','$mzPeaks','$intPeaks')";
	 echo($sql."\n");
	error_reporting(1);
	mysqli_query($link, $sql);
	$j++;
	
   }

   echo("Now hold on there pal....just about done...\nNow Moving to next library. \n");
   

 
   $i++;
}

 mysqli_close($conn);


?>
