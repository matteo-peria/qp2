Transcorrelated SCF
===================

### Three-electron contribution

With `three_body_h_tc=True`, each spin Fock matrix includes its
three-electron contribution. The total Fock matrix combines the spin
matrices and their three-electron terms are therefore already present.
With the default `three_body_h_tc=False` and `noL_standard=True`,
TC-SCF includes the scalar three-electron energy in its reported total
and omits the explicit three-electron term from the orbital update.
