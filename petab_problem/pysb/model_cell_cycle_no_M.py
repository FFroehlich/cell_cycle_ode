from pysb import (
    Model, Monomer, Parameter, Expression, Observable, Rule, Initial
)

model = Model('cell_cycle_no_M')

C = Monomer('live_cell',
            sites=['cycle'],
            site_states={'cycle': ['G1', 'G2', 'S']})

D = Monomer('dead_cell')

Observable('G1_obs', C(cycle='G1'))
Observable('S_obs', C(cycle='S'))
Observable('G2_plus_M_obs', C(cycle='G2'))
Observable('D_obs', D())

Rule('G1_S',    C(cycle='G1') >> C(cycle='S'),
     Expression('rateG1S',
                Parameter('kG1S', 0.1) * Parameter('rG1S', 0.1)))

Rule('S_G2',     C(cycle='S') >> C(cycle='G2'),
     Expression('rateSG2',
                Parameter('kSG2', 0.1) * Parameter('rSG2', 0.1)))

Rule('G2_M_G1', C(cycle='G2') >> C(cycle='G1') + C(cycle='G1'),
     Expression('rateG2MG1',
                Parameter('kG2MG1', 0.1) * Parameter('rG2MG1', 0.1)))

Rule('death', C() >> D(),
     Expression('ratedeath',
                Parameter('kphi', 0.1) * Parameter('rphi', 0.1)))

Initial(C(cycle='G1'), Parameter('G1_0'))
Initial(C(cycle='S'), Parameter('S_0'))
Initial(C(cycle='G2'), Parameter('G2_0'))
Initial(D(), Parameter('D_0'))
