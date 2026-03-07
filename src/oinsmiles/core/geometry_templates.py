"""
Geometry templates for OIN Engine.
Defines vector positions and reference vectors for standard geometries.
"""
import numpy as np
from typing import Dict, Any

# Structure: GEO -> Slot Idx -> {'pos': [x,y,z], 'ref': [x,y,z]}
# pos = The ideal slot vector (where the user puts the ligand)
# ref = The reference vector for North Star calculation (X_ref)

TEMPLATES: Dict[str, Dict[int, Dict[str, np.ndarray]]] = {
    'LIN': { # Linear
        0: {'pos': np.array([0,0,1]),  'ref': np.array([1,0,0])},
        1: {'pos': np.array([0,0,-1]), 'ref': np.array([1,0,0])}
    },
    'SPL': { # Square Planar (xy plane)
        0: {'pos': np.array([1,0,0]),  'ref': np.array([0,0,1])},
        1: {'pos': np.array([0,1,0]),  'ref': np.array([0,0,1])},
        2: {'pos': np.array([-1,0,0]), 'ref': np.array([0,0,1])},
        3: {'pos': np.array([0,-1,0]), 'ref': np.array([0,0,1])}
    },
    'OCT': { # Octahedral
        0: {'pos': np.array([0,0,1]),  'ref': np.array([1,0,0])}, # Axial Top
        1: {'pos': np.array([0,0,-1]), 'ref': np.array([1,0,0])}, # Axial Bottom
        2: {'pos': np.array([1,0,0]),  'ref': np.array([0,0,1])}, # Eq
        3: {'pos': np.array([-1,0,0]), 'ref': np.array([0,0,1])}, # Eq
        4: {'pos': np.array([0,1,0]),  'ref': np.array([0,0,1])}, # Eq
        5: {'pos': np.array([0,-1,0]), 'ref': np.array([0,0,1])}  # Eq
    },
    'TET': { # Tetrahedral
        0: {'pos': np.array([1,1,1]),   'ref': np.array([0,0,1])},
        1: {'pos': np.array([1,-1,-1]), 'ref': np.array([0,0,1])},
        2: {'pos': np.array([-1,1,-1]), 'ref': np.array([0,0,1])},
        3: {'pos': np.array([-1,-1,1]), 'ref': np.array([0,0,1])}
    },
    'TBP': { # Trigonal Bipyramidal
        0: {'pos': np.array([0,0,1]),            'ref': np.array([1,0,0])}, # Axial
        1: {'pos': np.array([0,0,-1]),           'ref': np.array([1,0,0])}, # Axial
        2: {'pos': np.array([0,1,0]),            'ref': np.array([0,0,1])}, # Eq
        3: {'pos': np.array([0.8660254,-0.5,0]), 'ref': np.array([0,0,1])}, # Eq
        4: {'pos': np.array([-0.8660254,-0.5,0]),'ref': np.array([0,0,1])}  # Eq
    },
    'TPY': { # Trigonal Pyramidal
        0: {'pos': np.array([0,0,1]),            'ref': np.array([1,0,0])},
        1: {'pos': np.array([0,1,0]),            'ref': np.array([0,0,1])},
        2: {'pos': np.array([0.8660254,-0.5,0]), 'ref': np.array([0,0,1])},
        3: {'pos': np.array([-0.8660254,-0.5,0]),'ref': np.array([0,0,1])}
    },
    'SPY': { # Square Pyramidal
        0: {'pos': np.array([0,0,1]),  'ref': np.array([1,0,0])}, # Axial
        1: {'pos': np.array([1,0,0]),  'ref': np.array([0,0,1])}, # Eq
        2: {'pos': np.array([-1,0,0]), 'ref': np.array([0,0,1])}, # Eq
        3: {'pos': np.array([0,1,0]),  'ref': np.array([0,0,1])}, # Eq
        4: {'pos': np.array([0,-1,0]), 'ref': np.array([0,0,1])}  # Eq
    },
    'TPL': { # Trigonal Planar
        0: {'pos': np.array([0,1,0]),            'ref': np.array([0,0,1])},
        1: {'pos': np.array([0.8660254,-0.5,0]), 'ref': np.array([0,0,1])},
        2: {'pos': np.array([-0.8660254,-0.5,0]),'ref': np.array([0,0,1])}
    },
    'PBP': { # Pentagonal Bipyramidal
        0: {'pos': np.array([0,0,1]), 'ref': np.array([1,0,0])}, # Axial
        1: {'pos': np.array([0,0,-1]), 'ref': np.array([1,0,0])}, # Axial
        2: {'pos': np.array([1,0,0]), 'ref': np.array([0,0,1])},
        3: {'pos': np.array([0.30901699, 0.95105652, 0]), 'ref': np.array([0,0,1])},
        4: {'pos': np.array([-0.80901699, 0.58778525, 0]), 'ref': np.array([0,0,1])},
        5: {'pos': np.array([-0.80901699, -0.58778525, 0]), 'ref': np.array([0,0,1])},
        6: {'pos': np.array([0.30901699, -0.95105652, 0]), 'ref': np.array([0,0,1])}
    }
}
