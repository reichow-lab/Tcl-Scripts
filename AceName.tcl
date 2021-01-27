proc g2a {inprefix} {
  mol new $inprefix.psf
  mol addfile $inprefix.pdb
  set all [atomselect top all]
  set chainatoms [atomselect top "protein and resid 2 and name CA"]
  set chainlist [$chainatoms get chain]
  foreach chain $chainlist {
    set ace [atomselect top "chain $chain and name CAY HY1 HY2 HY3 CY OY"]
    $ace set resid 1
    $ace set resname ACE
  }
  $all writepsf $inprefix.mod.psf
}
