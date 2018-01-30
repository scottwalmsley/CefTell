<?php

/**
 * class sampleMetaObj
 * A class for holding experiment metadata info for the experiment db.
 * Copyright 2014 Scott Walmsley, PhD
 *
 */


class sampleMetaObj{

	public $fileName;
	public $taxa;
	public $sampleType;
	public $polarity;
	public $sampleID;
	public $experimentID;
	public $samplePrep;
	public $extractionType;
	public $columnType;
	public $scanLevel;
	public $sampleRepID;
	public $technicalRepID;
	public $injectionRep;
	public $userDefinedRepID;
	#private $combinedCefFileName;

	public function setObjData($dataLine){
		$dataArr = explode("\t",$dataLine);
		$this->fileName  = $dataArr[0];
		$this->taxa = $dataArr[1];
		$this->sampleType = $dataArr[2];
		$this->polarity = $dataArr[3];
		$this->sampleID = $dataArr[4];
		$this->samplePrep = $dataArr[6];
		$this->extractionType = $dataArr[7];
		$this->columnType = $dataArr[8];
		$this->scanLevel = $dataArr[9];
		$this->sampleRepID = $dataArr[11];
		$this->technicalRepID = $dataArr[12];
		$this->injectionRep = $dataArr[10];
		#$this->combinedCefFileName = $dataArr[6];

	}
	public function getFileName() {
		return $this->fileName;
	}
	public function setFileName($fileName) {
		$this->fileName = $fileName;
		return $this;
	}
	public function getExperimentID() {
		return $this->experimentID;
	}
	public function setExperimentID($experimentID) {
		$this->experimentID = $experimentID;
		return $this;
	}
	public function getExtractionType() {
		return $this->extractionType;
	}
	public function setExtractionType($extractionType) {
		$this->extractionType = $extractionType;
		return $this;
	}
	public function getColumnType() {
		return $this->columnType;
	}
	public function setColumnType($columnType) {
		$this->columnType = $columnType;
		return $this;
	}
	public function getSampleRepID() {
		return $this->sampleRepID;
	}
	public function setSampleRepID($sampleRepID) {
		$this->sampleRepID = $sampleRepID;
		return $this;
	}
	public function getTechnicalRepID() {
		return $this->technicalRepID;
	}
	public function setTechnicalRepID($technicalRepID) {
		$this->technicalRepID = $technicalRepID;
		return $this;
	}
	public function getInjectionRep() {
		return $this->injectionRep;
	}
	public function setInjectionRep($injectionRep) {
		$this->injectionRep = $injectionRep;
		return $this;
	}
	public function getUserDefinedRepID() {
		return $this->userDefinedRepID;
	}
	public function setUserDefinedRepID($userDefinedRepID) {
		$this->userDefinedRepID = $userDefinedRepID;
		return $this;
	}
	public function getTaxa() {
		return $this->taxa;
	}
	public function setTaxa($taxa) {
		$this->taxa = $taxa;
		return $this;
	}
	public function getSampleType() {
		return $this->sampleType;
	}
	public function setSampleType($sampleType) {
		$this->sampleType = $sampleType;
		return $this;
	}
	public function getPolarity() {
		return $this->polarity;
	}
	public function setPolarity($polarity) {
		$this->polarity = $polarity;
		return $this;
	}
	public function getSampleID() {
		return $this->sampleID;
	}
	public function setSampleID($sampleID) {
		$this->sampleID = $sampleID;
		return $this;
	}
	public function getSamplePrep() {
		return $this->samplePrep;
	}
	public function setSamplePrep($samplePrep) {
		$this->samplePrep = $samplePrep;
		return $this;
	}
	public function getScanLevel() {
		return $this->scanLevel;
	}
	public function setScanLevel($scanLevel) {
		$this->scanLevel = $scanLevel;
		return $this;
	}




}


?>
