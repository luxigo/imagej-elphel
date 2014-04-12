<?php
/*!*******************************************************************************
*! FILE NAME  : escher_pattern.php
*! DESCRIPTION: slanted curved checker board generator
*! requires php5-ps, ps2pdf
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
$tmpdir ='/tmp';
$filetype='PDF';

foreach ($_GET as $key=>$value) {
  switch (strtoupper($key)) {
    case 'PAGE_WIDTH': $PAGE_WIDTH=$value+0;break;
    case 'PAGE_HEIGHT':$PAGE_HEIGHT=$value+0;break;
    case 'LPM':        $LPM=$value+0;break;
    case 'ESCHER':     $ESCHER=$value+0;break;
    case 'ROTATE':     $ROTATE=$value+0;break;
    case 'MODE':
    case 'TYPE':       $filetype=$value;break;
  }
}

$basename='escher-pattern-ESCHER'.$ESCHER.'-LPM'.$LPM.'-ROT'.$ROTATE.'-PAGE_WIDTH'.$PAGE_WIDTH.'-PAGE_HEIGHT'.$PAGE_HEIGHT;
$ps_name=$tmpdir.'/'.$basename.'.ps';
$pdf_name=$tmpdir.'/'.$basename.'.pdf';

$mm2points=72/25.4;
$PW_POINTS=$PAGE_WIDTH*$mm2points;
$PH_POINTS=$PAGE_HEIGHT*$mm2points;

   $ps = ps_new();
   if (!ps_open_file($ps, $ps_name)) {
     print 'Cannot open PostScript file '.$ps_name."\n";
     exit(1);
   }
   ps_set_parameter($ps, "warning", "true");
   ps_set_info($ps, "Creator", " escher_pattern.php");
   ps_set_info($ps, "Author", "Elphel");
   ps_set_info($ps, "Title", "PSF measurement pattern");

   $Size=1.0;
   $qSize=$Size/4;
   $hSize=$Size/2;
   $a=$ESCHER*(sqrt(2)-1.0);
   $angle=-$ROTATE; // PS starts from tyhe bottom left
   if ($ESCHER>0) {
     $r=($a*$a+1)/(2*$a)*$qSize;
     $r2=$r*$r;
     $h=sqrt($r2-$qSize*$qSize);
     $dc=2*$qSize-$h;
     $halfAangle=rad2deg(atan($qSize/$h));
//$center=1.5*$sideOfSquare;
     $center=$dc+$r;
     $pstemplate = ps_begin_template($ps, 2*$center, 2*$center);
       ps_moveto($ps, $center-$hSize, $center-$hSize);
       ps_arc   ($ps  , $center-$Size+$dc,   $center-$qSize,   $r, 0-$halfAangle  ,   0+$halfAangle  );
       ps_arcn  ($ps  , $center-$dc  ,    $center+$qSize,    $r, 180+$halfAangle  , 180-$halfAangle  );
       ps_arc   ($ps  , $center-$qSize,   $center+$Size-$dc, $r, 270-$halfAangle  , 270+$halfAangle  );
       ps_arcn  ($ps  , $center+$qSize,   $center+$dc,       $r, 90+$halfAangle  ,   90-$halfAangle  );
       ps_arc   ($ps  , $center+$Size-$dc,$center+$qSize,    $r, 180-$halfAangle  , 180+$halfAangle  );
       ps_arcn  ($ps  , $center+$dc,    $center-$qSize,      $r,   0+$halfAangle  ,   0-$halfAangle  );
       ps_arc   ($ps  , $center+$qSize,   $center-$Size+$dc, $r,  90-$halfAangle  ,  90+$halfAangle  );
       ps_arcn  ($ps  , $center-$qSize,   $center-$dc,       $r, 270+$halfAangle  , 270-$halfAangle  );
       ps_fill($ps);
     ps_end_template($ps);
   } else {
     $pstemplate = ps_begin_template($ps, $Size, $Size);
       ps_moveto($ps, 0,    0);
       ps_lineto($ps, 0,    $Size);
       ps_lineto($ps, $Size,$Size);
       ps_lineto($ps, $Size,0);
       ps_lineto($ps, 0,0);
       ps_fill($ps);
     ps_end_template($ps);
   }

   $s=$PW_POINTS+$PH_POINTS;
   $sideOfSquare=500/$LPM*$mm2points;
   $numberOfSquares=floor((4*$s-$sideOfSquare)/(2*$sideOfSquare));
   $s=(2*$numberOfSquares+1)*$sideOfSquare;

   ps_begin_page($ps, $PW_POINTS, $PH_POINTS);

   ps_setcolor($ps, "fill", "cmyk", 0, 0, 0, 1);
   ps_rotate($ps,$angle);


   for ($y=-$s;$y<$s;$y+=2*$sideOfSquare) for ($x=-$s;$x<$s;$x+=2*$sideOfSquare) {
     ps_place_image($ps, $pstemplate, $x, $y,$sideOfSquare);
     ps_place_image($ps, $pstemplate, $x+$sideOfSquare, $y+$sideOfSquare,$sideOfSquare);
   }
   ps_end_page($ps);
   ps_close($ps);
   ps_delete($ps);
   if ($filetype=='PDF') {
//   echo 'ps2pdf -dDEVICEWIDTHPOINTS='.$PW_POINTS.' -dDEVICEHEIGHTPOINTS='.$PH_POINTS.' '.$ps_name.' '.$pdf_name;
     exec('ps2pdf -dDEVICEWIDTHPOINTS='.$PW_POINTS.' -dDEVICEHEIGHTPOINTS='.$PH_POINTS.' '.$ps_name.' '.$pdf_name);
     $outfile=$pdf_name;
     header('Content-Type: application/pdf');
  } else {
     $outfile=$ps_name;
     header('Content-Type: application/postscript');
  }
   if (headers_sent()) {
     echo 'Some data has already been output to browser, can\'t send PDF file';
     exit (1);
   }
   header('Cache-Control: public, must-revalidate, max-age=0'); // HTTP/1.1
   header('Pragma: public');
   header('Expires: Sat, 26 Jul 1997 05:00:00 GMT'); // Date in the past
   header('Last-Modified: '.gmdate('D, d M Y H:i:s').' GMT');
   $fp = fopen($outfile, 'rb');
   fseek($fp, 0, SEEK_END);  /// file pointer at the end of the file (to find the file size)
   $fsize = ftell($fp);      /// get file size
   fseek($fp, 0, SEEK_SET);  /// rewind to the start of the file
   header('Content-Length: '.$fsize);
   header('Content-Disposition: inline; filename="'.basename($outfile).'";');
   fpassthru($fp);           /// send the raw data itself
   fclose($fp);

//echo 'Done - see '.$pdf_name;

                
?> 
