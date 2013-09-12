<?php include "head.php" ; ?>
<html>
<body>
<?php include "topmenu.php" ; ?>
<?php include "imgbar.php" ; ?>

<div id="page">
   <div id="content_full">
      <div class="post hr">

         <h2>Comet parameter: peptide_mass_units</h2>

         <ul>
         <li>This parameter controls the units applied to the <a href="peptide_mass_tolerance.php">peptide_mass_tolerance</a> parameter.
         <li>Valid values are 0, 1, and 2:
            <ul>
            <li>0 for amu
            <li>1 for mmu
            <li>2 for ppm
            </ul>
         </ul>

         <p>Example:
         <br><tt>peptide_mass_units = 0</tt>
         <br><tt>peptide_mass_units = 1</tt>
         <br><tt>peptide_mass_units = 2</tt>

      </div>
   </div>
   <div style="clear: both;">&nbsp;</div>
</div>

<?php include "footer.php" ; ?>
</html>
