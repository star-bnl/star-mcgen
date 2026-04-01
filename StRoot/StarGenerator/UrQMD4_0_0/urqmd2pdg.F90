integer function urqmd_pdgid( ityp, iso3 ) bind (c, name="urqmd_pdgid_" )
  integer, intent(in) :: ityp, iso3
  integer, external :: pdgid
  urqmd_pdgid = pdgid(ityp,iso3)
end function urqmd_pdgid
