<?php
/*!*******************************************************************************
*! FILE NAME  : slanted_checker.php
*! DESCRIPTION: slanted checker board target
*! Copyright (C) 2010 Elphel, Inc
*! -----------------------------------------------------------------------------**
*!
*!  This program is free software: you can redistribute it and/or modify
*!  it under the terms of the GNU General Public License as published by
*!  the Free Software Foundation, either version 3 of the License, or
*!  (at your option) any later version.
*!
*!  This program is distributed in the hope that it will be useful,
*!  but WITHOUT ANY WARRANTY; without even the implied warranty of
*!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*!  GNU General Public License for more details.
*!
*!  You should have received a copy of the GNU General Public License
*!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*! -----------------------------------------------------------------------------**
*/

require_once('../tcpdf/config/lang/eng.php');
require_once('../tcpdf/tcpdf.php');
$orientation=   'P'; // P/L
$unit=          'mm' ; /// pt,mm,cm,in
//$format=        array(40,20); /// width and height of current units
$unicode=       true;
$encoding=      "UTF-8";
$diskcache=     false; /// use temporary files
$margin_top=    2;
$margin_bottom= 1;
$margin_left=   1;
$margin_right=  1;
$LPM=50;
$PAGE_WIDTH=270;
$PAGE_HEIGHT=210;
$ROTATE=5;
$ESCHER=2.0;
$NC=3; // number of curves for the circle (seems 3 is OK for small images, 2 - visible difference)
foreach ($_GET as $key=>$value) {
  switch (strtoupper($key)) {
    case 'PAGE_WIDTH': $PAGE_WIDTH=$value+0;break;
    case 'PAGE_HEIGHT':$PAGE_HEIGHT=$value+0;break;
    case 'LPM':        $LPM=$value+0;break;
    case 'ESCHER':     $ESCHER=$value+0;break;
    case 'NC':         $NC=$value+0;break;
    case 'ROTATE':     $ROTATE=$value+0;break;
  }
}

$format=        array($PAGE_WIDTH,$PAGE_HEIGHT); /// width and height of current units 279.4 x 215.9 (11.5x8.5)




$pdf=new TCPDF ($orientation,
                $unit,
                $format,
                $unicode,
                $encoding,
                $diskcache);

// set document information
$pdf->SetCreator(PDF_CREATOR);
$pdf->SetAuthor('Elphel');
$pdf->SetTitle('PSF chart');


/// set default header data
/// remove default header/footer
$pdf->setPrintHeader(false);
$pdf->setPrintFooter(false);

/// set default monospaced font
$pdf->SetDefaultMonospacedFont(PDF_FONT_MONOSPACED);

///set margins
$pdf->SetMargins($margin_left, $margin_top, $margin_right);
//set auto page breaks
$pdf->SetAutoPageBreak(FALSE, $margin_bottom);
$pdf->AddPage();


   $pdf->StartTransform();
   $h=$format[1]/2;
   $w=$format[0]/2;
   $pdf->Rotate(-$ROTATE,$w,$h);
   //! $ROTATE>0!

   $s=$w+$h;
   $sideOfSquare=500/$LPM;
   $numberOfSquares=floor((4*$s-$sideOfSquare)/(2*$sideOfSquare));
   $s=(2*$numberOfSquares+1)*$sideOfSquare;
   for ($y=-$s;$y<$s;$y+=2*$sideOfSquare) for ($x=-$s;$x<$s;$x+=2*$sideOfSquare) {
        $pdf->Rect($x,               $y,               $sideOfSquare, $sideOfSquare,'F',array(),array(0,0,0)); 
        $pdf->Rect($x+$sideOfSquare, $y+$sideOfSquare, $sideOfSquare, $sideOfSquare,'F',array(),array(0,0,0)); 
   }

   if ($ESCHER>0) {
      $qSize=$sideOfSquare/4;
      $a=$ESCHER*(sqrt(2)-1.0);
      $r=($a*$a+1)/(2*$a)*$qSize;
      $r2=$r*$r;
      $h=sqrt($r2-$qSize*$qSize);
      $dc=2*$qSize-$h;
      $halfAangle=rad2deg(atan($qSize/$h));
      for ($y=-$s;$y<$s;$y+=2*$sideOfSquare) for ($x=-$s;$x<$s;$x+=2*$sideOfSquare) {
        for ($subY=0;$subY<2;$subY++) for ($subX=0;$subX<2;$subX++){
          $cellColor=($subY==$subX)?array(0,0,0):array(255,255,255);
          $xc=$x+$sideOfSquare*$subX+2*$qSize;
          $yc=$y+$sideOfSquare*$subY+2*$qSize;
/** Limited precision make nasty lines between circle segments , so now they are always semi-circles) */
          $pdf->Circle($xc-$dc,    $yc-$qSize, $r, 180-$halfAangle,    -$halfAangle, 'F', array(), $cellColor,$NC); // , $nc=8 - number of curves - last parameter
          $pdf->Circle($xc+$qSize, $yc-$dc,    $r,  90-$halfAangle, 270-$halfAangle, 'F', array(), $cellColor,$NC); // , $nc=8 - number of curves - last parameter
          $pdf->Circle($xc+$dc,    $yc+$qSize, $r,    -$halfAangle, 180-$halfAangle, 'F', array(), $cellColor,$NC); // , $nc=8 - number of curves - last parameter
          $pdf->Circle($xc-$qSize, $yc+$dc,    $r, 270-$halfAangle,  90-$halfAangle, 'F', array(), $cellColor,$NC); // , $nc=8 - number of curves - last parameter

        }
      }
   }


   $pdf->StopTransform();


$pdf->Output('slanted-'.(($ESCHER>0)?('ESCHER-'.$ESCHER):'').'-'.$LPM.'LPM-'.$ROTATE.'grad.pdf', 'I');
                
?> 
