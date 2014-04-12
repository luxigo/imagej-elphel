<?php
/*!*******************************************************************************
*! FILE NAME  : target_circle.php
*! DESCRIPTION: Creates target of concentric circles
*! REQUIRES: tcpdf
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

require_once('tcpdf/config/lang/eng.php'); /// external dependency
require_once('tcpdf/tcpdf.php'); /// external dependency
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
$TARGET_DIAMETER=100;
$TARGET_CIRCLES= 5;
$CENTER_WHITE=   1.0  ; //put white dot in the center

foreach ($_GET as $key=>$value) {
  switch (strtoupper($key)) {
    case 'D':  
    case 'TARGET_DIAMETER':   $TARGET_DIAMETER=$value+0;break;
    case 'C':  
    case 'N':  
    case 'TARGET_CIRCLES':    $TARGET_CIRCLES=$value+0;break;
    case 'W':  
    case 'CENTER_WHITE':    $CENTER_WHITE=$value+0;break;

  }
}


 $PAGE_WIDTH= ceil($TARGET_DIAMETER+ $margin_left+$margin_right);
 $PAGE_HEIGHT=ceil($TARGET_DIAMETER+$margin_top+ $margin_bottom);

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
 $pdf->SetTitle('Focusing Target');


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
 for ($i=$TARGET_CIRCLES*2; $i>=0; $i--) {
   $r=($i+1)*$TARGET_DIAMETER/($TARGET_CIRCLES*2+1)/2;
   
   $pdf->Circle($TARGET_DIAMETER/2,$TARGET_DIAMETER/2,$r,0,360,'FD', array(), ($i & 1)? array(255,255,255):array(BLACK));
 }
 if ($CENTER_WHITE>0) {
   $pdf->Circle($TARGET_DIAMETER/2,$TARGET_DIAMETER/2,$CENTER_WHITE/2,0,360,'FD', array(),  array(255,255,255));
 }

 $pdf->Output('target_circles_D'.$TARGET_DIAMETER.'_N'.$TARGET_CIRCLES.'.pdf', 'I');
 exit(0);

		/**
		* Draws a circle.
		* A circle is formed from n Bezier curves.
		* @param float $x0 Abscissa of center point.
		* @param float $y0 Ordinate of center point.
		* @param float $r Radius.
		* @param float $astart: Angle start of draw line. Default value: 0.
		* @param float $afinish: Angle finish of draw line. Default value: 360.
		* @param string $style Style of rendering. Possible values are:
		* <ul>
		*	 <li>D or empty string: Draw (default).</li>
		*	 <li>F: Fill.</li>
		*	 <li>DF or FD: Draw and fill.</li>
		*	 <li>C: Draw close.</li>
		*	 <li>CNZ: Clipping mode (using the even-odd rule to determine which regions lie inside the clipping path).</li>
		*	 <li>CEO: Clipping mode (using the nonzero winding number rule to determine which regions lie inside the clipping path).</li>
		* </ul>
		* @param array $line_style Line style of circle. Array like for {@link SetLineStyle SetLineStyle}. Default value: default line style (empty array).
		* @param array $fill_color Fill color. Format: array(red, green, blue). Default value: default color (empty array).
		* @param integer $nc Number of curves used in circle. Default value: 8.
		* @access public
		* @since 2.1.000 (2008-01-08)
		*/
///		public function Circle($x0, $y0, $r, $astart=0, $afinish=360, $style='', $line_style=array(), $fill_color=array(), $nc=8) {

                
?> 
