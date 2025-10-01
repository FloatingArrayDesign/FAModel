import numpy as np
import moorpy as mp
from famodel.design.LineDesign import LineDesign
from pathlib import Path



moorpy_dir = Path(mp.__file__).parent
moorprops_file = moorpy_dir / "library" / "MoorProps_default.yaml"
moorprops_file = str(moorprops_file)


def test_DEA_Chain_COBYLA_none():

    depth = 200

    settings = {}
    settings['rBFair'] = [58,0,-14]
    settings['x_ampl'] = 10
    settings['fx_target'] = 1.95e6 + 465284.3 # from 0.71^2 m/s * 9.23e5 current
    settings['kx_target'] = 2e4
    settings['headings'] = [60, 180, 300]

    settings['name'] = 'DEA-chain'
    settings['lineTypeNames'] = ['chain']
    settings['anchorType'] = 'drag-embedment'
    settings['allVars'] = [1000/10, 1000, 120]

    # solve_for = none setup
    settings['solve_for'] = 'none'
    settings['x_target'] = 35
    settings['x_mean_out'] = 35
    settings['x_mean_in'] = 60

    settings['Xindices'] = [0, 1, 2]
    settings['Xmin'] = [10, 10, 10]
    settings['Xmax'] = [1000, 1000, 120]
    settings['dX_last'] = [10, 10, 10]

    settings['constraints'] = [dict(name='min_lay_length', index=0, threshold=20, offset='max'),
                                dict(name='max_offset'    , index=0, threshold=60, offset='min'),
                                dict(name='max_offset'    , index=0, threshold=35, offset='max')]
    for j in range(len(settings['lineTypeNames'])):
        settings['constraints'].append(dict(name='tension_safety_factor', index=j, threshold=2.0, offset='max'))


    # run LineDesign
    ld1 = LineDesign(depth, lineProps=moorprops_file, **settings)
    ld1.setNormalization()
    X, min_cost = ld1.optimize(maxIter=400, plot=False, display=2, stepfac=4, method='COBYLA')
    ld1.updateDesign(X, display=0)
    X1 = np.array(X)*ld1.X_denorm
    X1[0] *= 10

    assert abs(X1[0] - 1033.070) < 0.01
    assert abs(X1[1] - 1036.494) < 0.01
    assert abs(X1[2] -  102.920) < 0.01



def test_DEA_Chain_COBYLA_offset():

    depth = 200

    settings = {}
    settings['rBFair'] = [58,0,-14]
    settings['x_ampl'] = 10
    settings['fx_target'] = 1.95e6 + 465284.3 # from 0.71^2 m/s * 9.23e5 current
    settings['kx_target'] = 2e4
    settings['headings'] = [60, 180, 300]

    settings['name'] = 'DEA-chain'
    settings['lineTypeNames'] = ['chain']
    settings['anchorType'] = 'drag-embedment'
    settings['allVars'] = [1000/10, 1000, 120]

    settings['solve_for'] = 'offset'
    settings['x_target'] = 34.4
    settings['x_mean_out'] = 34.4
    settings['x_mean_in'] = 60

    settings['Xindices'] = [0, 's', 1]
    settings['Xmin'] = [10, 10]
    settings['Xmax'] = [500, 500]
    settings['dX_last'] = [10, 10]

    settings['constraints'] = [dict(name='min_lay_length', index=0, threshold=20, offset='max')]
    for j in range(len(settings['lineTypeNames'])):
        settings['constraints'].append(dict(name='tension_safety_factor', index=j, threshold=2.0, offset='max'))


    # run LineDesign
    ld2 = LineDesign(depth, lineProps=moorprops_file, **settings)
    ld2.setNormalization()
    X, min_cost = ld2.optimize(maxIter=400, plot=False, display=2, stepfac=4, method='COBYLA')
    ld2.updateDesign(X, display=0)
    X2 = np.array(X)*ld2.X_denorm
    X2[0] *= 10

    assert abs(X2[0] - 1033.144) < 0.01
    assert abs(ld2.ss.lineList[0].L - 1036.553) < 0.01
    assert abs(X2[1] -  102.869) < 0.01


def test_DEA_Chain_COBYLA_tension():

    depth = 200

    settings = {}
    settings['rBFair'] = [58,0,-14]
    settings['x_ampl'] = 10
    settings['fx_target'] = 1.95e6 + 465284.3 # from 0.71^2 m/s * 9.23e5 current
    settings['kx_target'] = 2e4
    settings['headings'] = [60, 180, 300]

    settings['name'] = 'DEA-chain'
    settings['lineTypeNames'] = ['chain']
    settings['anchorType'] = 'drag-embedment'
    settings['allVars'] = [1000/10, 1000, 120]

    settings['solve_for'] = 'tension'
    settings['x_target'] = 34.4
    settings['x_mean_out'] = 34.4
    settings['x_mean_in'] = 60
    settings['fx_target'] = 570217.4

    settings['Xindices'] = [0, 's', 1]
    settings['Xmin'] = [10, 10]
    settings['Xmax'] = [500, 500]
    settings['dX_last'] = [10, 10]

    settings['constraints'] = [dict(name='min_lay_length', index=0, threshold=20, offset='max')]
    for j in range(len(settings['lineTypeNames'])):
        settings['constraints'].append(dict(name='tension_safety_factor', index=j, threshold=2.0, offset='max'))


    # run LineDesign
    ld3 = LineDesign(depth, lineProps=moorprops_file, **settings)
    ld3.setNormalization()
    X, min_cost = ld3.optimize(maxIter=400, plot=False, display=2, stepfac=4, method='COBYLA')
    ld3.updateDesign(X, display=0)
    X3 = np.array(X)*ld3.X_denorm
    X3[0] *= 10

    assert abs(X3[0] - 1033.145) < 0.01
    assert abs(ld3.ss.lineList[0].L - 1036.553) < 0.01
    assert abs(X3[1] -  102.869) < 0.01



'''
def test_DEA_Chain_Polyester_COBYLA_none():

    depth = 200

    settings = {}
    settings['rBFair'] = [58,0,-14]
    settings['x_ampl'] = 10
    settings['fx_target'] = 1.95e6 + 465284.3 # from 0.71^2 m/s * 9.23e5 current
    settings['headings'] = [60, 180, 300]

    settings['solve_for'] = 'none'
    settings['x_target'] = 35
    settings['x_mean_out'] = 35
    settings['x_mean_in'] = 60

    settings['name'] = 'DEA-chain-rope'
    settings['lineTypeNames'] = ['chain','polyester']
    settings['anchorType'] = 'drag-embedment'

    settings['solve_for'] = 'none'
    settings['allVars'] = [750/10, 600, 100, 0, 100, 180]
    settings['Xindices'] = [0, 1, 2, 'c', 3, 4]
    settings['Xmin'] = [10, 10, 10, 10, 10]
    settings['Xmax'] = [500, 1500, 300, 1000, 300]
    settings['dX_last'] = [10, 10, 10, 10, 10]

    settings['constraints'] = [dict(name='min_lay_length', index=0, threshold=20, offset='max'),
                            dict(name='max_offset'    , index=0, threshold=60, offset='min'),
                            dict(name='max_offset'    , index=0, threshold=35, offset='max'),
                            dict(name='rope_contact'  , index=1, threshold=5 , offset='min')]
    for j in range(len(settings['lineTypeNames'])):
        settings['constraints'].append(dict(name='tension_safety_factor', index=j, threshold=2.0, offset='max'))

    ld = LineDesign(depth, lineProps=moorprops_file, **settings)
    ld.setNormalization()
    X, min_cost = ld.optimize(maxIter=4000, plot=False, display=3, stepfac=4, method='COBYLA')
    ld.updateDesign(X, display=0)
    X = np.array(X)*ld1.X_denorm
    X[0] *= 10

    assert abs(X[0] - 986.242) < 0.001
    assert abs(X[1] - 856.697) < 0.001
    assert abs(X[2] - 106.746) < 0.001
    assert abs(X[3] - 131.426) < 0.001
    assert abs(X[4] - 179.636) < 0.001
'''



