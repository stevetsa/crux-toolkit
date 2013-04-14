<?php include "head.php" ; ?>
<html>
<body>
<?php include "topmenu.php" ; ?>
<?php include "imgbar.php" ; ?>

<div id="page">
   <div id="content_full">
      <div class="post hr">

         <h2>Comet parameter: minimum_peaks</h2>

         <ul>
         <li>An integer value indicating the minimum number of m/z-intensity pairs
         that must be present in a spectrum before it is searched.
         <li>This parameter can be used to avoid searching nearly sparse spectra
         that will not likely yield an indentification.
         <li>Valid values are any integer number.
         </ul>

         <p>Example:
         <br><tt>minimum_peaks = 20</tt>

      </div>
   </div>
   <div style="clear: both;">&nbsp;</div>
</div>

<?php include "footer.php" ; ?>
</html>
