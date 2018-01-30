<?php
/**
 * Class experimentDbMetaSQL
 * Class for building SQL tables for experiment evidence
 *Copyright 2014, Scott Walmsley, PhD
 *
 */
class experimentDbMetaSQL{
	
	/**
	 * 
	 * @param string $link
	 * @param string $tableName
	 */
	public static function createExptMetaSQLdb($link, $tableName){
		mysqli_query($link, 
			"CREATE TABLE $tableName (".
			"experimentID TINYTEXT,".
			"numberOfSamples INT(3),".
			"numberOfReplicates INT(3),".
			"numberOfInjections INT(3),".
			"extractionMethod TINYTEXT,".
			"instrumentMethod TINYTEXT,".
			"columnType TINYTEXT,".
			"polarity VARCHAR(8))"
		);
	}
	
	
	
	/**
	 * 
	 * @param string $link
	 * @param string $tableName
	 */
	public static function createSampleMetaSQLdb($link, $tableName){
		mysqli_query($link,
			"CREATE TABLE $tableName (".
			"sampleID VARCHAR(50),".
			"experimentID INT(5),".
			"technicalRep INT(3),".
			"sampleRep INT(3),".
			"injRep INT(3),".
			"extraction TINYTEXT,".
			"polarity TINYTEXT,".
			"col TINYTEXT,".
			"scanLevel INT(2),".
			"samplePrep TINYTEXT,".		
			"fileName VARCHAR(200))"
		);
	}
	
	/**
	 * 
	 * @param string $link
	 * @param string $tableName
	 */
	public static function createMSIonDataSQLdb($link,$tableName){
		mysqli_query($link,
			"CREATE TABLE $tableName (".
			"ionID INT(8) NOT NULL PRIMARYKEY,".
			"fileName TINYTEXT,".
			"compoundName TEXT,".
			"matchScore DOUBLE(7,3),".
			"molecularFormula TINYTEXT,".
			"precursorMZ DOUBLE(12,6),".
			"ppmError DOUBLE(7,3),".
			"retentionTime DOUBLE (10,3),".
			"qmsRT DOUBLE (4,3),".
			"scanIndex INT(9),".
			"hasMSMS VARCHAR(1),".
			"numberCEvs INT(2),".
			"numberMSMSevents INT(3),".
			"meanMSMSscore DOUBLE(4,1))"
		);
		
	}
	
	/**
	 *
	 * @param string $link
	 * @param string $tableName
	 */
	public static function createMSMSIonDataSQLdb($link,$tableName){
		mysqli_query($link,
			"CREATE TABLE $tableName (".
			"spectrumID INT(8),".
			"isUnknown VARCHAR(1),".
			"library TINYTEXT,".
			"libraryID VARCHAR(20),".
			"ionID INT(8),".
			"fileName TINYTEXT,".
			"isUnknown VARCHAR (1)".
			"dot INT (3)".
			"revDot INT(3)".
			"probability DOUBLE(4,1))"
		);	
	}
	
	
	
	
	
	
}

?>
