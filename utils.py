


import numpy as np
import quaternion

import os, subprocess


################################ Utility Methods ###############################

this_path = os.path.dirname(os.path.abspath(__file__))

def gen_PA(name):
    # Write a gen_PA.pgn file specific to PA_name

    PA_name = '%s_specific'%name

    data = [
    'package require psfgen',
    'psfcontext reset',
    'topology %s/top_all36_prot_forPAs.rtf'%this_path,
    'pdbalias residue HIS HSD',
    'segment P001 {',
    'pdb %s.pdb'%PA_name,
    'first none',
    'last none',
    'auto angles dihedrals',
    '}',
    'coordpdb %s.pdb P001'%PA_name,
    'regenerate angles dihedrals',
    'guesscoord',
    'writepdb %s_aa.pdb'%name,
    'exit'
    ]

    f = open('gen_%s.pgn'%name, 'wt')
    f.write('\n'.join(data))
    f.close


def generic_to_specific_PA(PA_seq, PA_generic, name):
    """Write up a specific PA sequence PA_seq : <name>_specific.pdb
    PA_seq is in string format. Example: C12VVAAEE. 
    Digits are used only for the alkyl chain.
    Use the amount of alkyl chain asked for. Currently the maximum is C16.
    Exchange A,B,C... in PA_generic.pdb to valid residues
    """
    
    if (PA_seq[0] == 'C') and PA_seq[1].isdigit():
        if PA_seq[2].isdigit():
            num_alkylC = int(PA_seq[1:3])
            pep_seq = PA_seq[3:]
        else:
            num_alkylC = int(PA_seq[1])
            pep_seq = PA_seq[2:]
    else: # no alkyl
        pep_seq = PA_seq
        num_alkylC = 0

    num_res = len(pep_seq)


    res_dict = {'A':"ALA",'C':"CYS",'D':"ASP",'E':"GLU",'F':"PHE", 
                'G':"GLY",'H':"HIS",'I':"ILE",'K':"LYS",'L':"LEU", 
                'M':"MET",'N':"ASN",'P':"PRO",'Q':"GLN",'R':"ARG", 
                'S':"SER",'T':"THR",'V':"VAL",'W':"TRP",'Y':"TYR"}

    
    res_dict['Z'] = '12C'
    res_dict['X'] = 'C12'


    # generic residue name in PA_generic.pdb
    generic_residues = ['AAA','BBB','CCC','DDD','EEE','FFF','GGG','HHH','III','JJJ',
                        'KKK','LLL','MMM','NNN','OOO','PPP','QQQ','RRR','SSS','TTT',
                        'UUU','VVV','WWW','XXX','YYY','ZZZ'] + \
                       ['AZZ','BZZ','CZZ','DZZ','EZZ','FZZ','GZZ','HZZ','IZZ','JZZ',
                        'KZZ','LZZ','MZZ','NZZ','OZZ','PZZ','QZZ','RZZ','SZZ','TZZ',
                       ] 
               
    
    with open(f'{this_path}/{PA_generic}', 'r') as f:
        lines = f.readlines()

    
    data = ''
    for line in lines:
        if ('C16' in line) and PA_seq[:3]!='C16':
            continue
        
        if generic_residues[num_res] in line:
            data = data +'END'
            break
        
        for i in range(num_res):
            if generic_residues[i] in line:
                line = line.replace(generic_residues[i], res_dict[pep_seq[i]])
                break

        data = data + line

    with open('%s_specific.pdb'%name, 'w') as fout:
        fout.write(data)


def make_aa_pdb(name):
    """Make the all-atom pdb file <name>_aa.pdb for the PA_seq using vmd
    
    """
    #generic_to_specific_PA(PA_seq.upper(), name)
    gen_PA(name)
    process = subprocess.run('vmd -dispdev text -e gen_%s.pgn'%name, shell=True)
    process.check_returncode()
    
    


def merge_topfiles(topfiles, newtopfile):
    """Merge multiple top files into one.
    Essentially appends the topfiles however
    - excluding #includes already included
    - changing the name of the protein
    """
    
    liness=[]
    for topfile in topfiles:
        with open(topfile, 'r') as f:
            liness += [f.readlines()]

    

    lines = []
    data = ''
    for i,lines_ in enumerate(liness):
        for line in lines_:
            if ('#include' in line) and (line in data):
                continue
            else:
                if 'Protein ' in line:
                    lines += [ line.replace('Protein ', 'Protein%s '%(i+1)) ]
                else:
                    lines += [ line ]
                data += line


    with open(newtopfile, 'w') as f:
        f.write(''.join(lines))



def update_num_molecules_in_topfile(molname, topfile, num_molecules):
    """
    Read the number of molecules added in the <grofile>
    and update that number in <topfile>
    """

    with open(topfile, 'r') as f:
        lines = f.readlines()
        start = False
        for i,line in enumerate(lines):
            if '#mols' in line:
                if lines[i+1].split()[0]==molname:
                    start = True
                    continue
            if start:
                lines[i] = f'{molname}       {num_molecules}\n'
                break
    
    with open(topfile, 'w') as f:
        f.write(''.join(lines))





def init_fiber_config(grofile, num_atoms, num_molecules, Lx, Ly, Lz, start_from_nth_atom=0, invert=False, C_indices=[0,2],
    delta = 0.5, radial_offset = 0.5):
    """ Change the random positioning in grofile with fiber like initial config
    assumes all molecules are in contiguous in the grofile
    Note: start_from_nth_atom starts from index 0, whereas in grofile atom index starts from 1
    Assumes first residue is C
    C_indices: Carbon indices to use when aligning. 
    Vector C_indices[1]-C_indices[0] is used.
    If invert is true, the other end is pointed towards fiber inside
    
    delta: distance between layers in nm
    radial_offset: distance form the center to inner PA tip 
    """


    # Parameters for the initial configuration as a fiber 
    num_PA_layer = 9
    theta = np.pi*40/180 # angle between two PA in a layer
    theta_offset = np.pi*20/180


    # Molecule's atoms' positions. Read from grofile
    atoms_positions = []
    with open(grofile, 'r') as f:
        lines = f.readlines()
        start = False
        C_pos = []
        for line in lines[2:]:
            words = line.split()
            try:
                float(words[-1])
            except:
                continue
            if len(words)<4:
                continue

            index = int(line[-30:-25])
            if index == start_from_nth_atom+1:
                start = True
            if start:
                if len(atoms_positions) == num_atoms:
                    break
                atoms_positions += [ [float(words[-3]), float(words[-2]), float(words[-1])] ]
                if 'C' == words[1][0]:
                    C_pos += [ atoms_positions[-1] ]
    
    if invert:
        atoms_positions = np.array(atoms_positions)[::-1]
    else:
        atoms_positions = np.array(atoms_positions)

    # Translate with origin at outermost (C16) Carbon of alkyl tail
    atoms_positions -= atoms_positions[C_indices[0]]

    # Align atoms_positions in the x-y plane by aligning the vector 
    # between first two consecutive C atoms of C16 towards x-axis
    # NOTE: Assumes the first atoms is C16
    v = np.array(C_pos[C_indices[1]]) - np.array(C_pos[C_indices[0]])
    q = quaternion.q_between_vectors(v, [1,0,0])
    for i,pos in enumerate(atoms_positions):
        atoms_positions[i] = quaternion.qv_mult(q, pos)
    
        
    C_positions = []  # of the outermost C
    quaternions = [] # wrt when PAM chain is oriented in x-axis
    for i in range(num_molecules):
        j = i  % num_PA_layer
        k = i // num_PA_layer
        
        q = quaternion.axisangle_to_q([0,0,1], theta * j + theta_offset * (1+(-1)**k)/2 )
        quaternions += [q]

        pos = np.array([0,0,k*delta]) + quaternion.qv_mult(q, [radial_offset,0,0])
        C_positions += [ pos ]


    # Calculate positions of all atoms
    positions = []
    for i in range(num_molecules):
        for j,pos in enumerate(atoms_positions):
            p = quaternion.qv_mult(quaternions[i], pos)
            p += C_positions[i]
            positions += [ p ]
    positions = np.array(positions) - np.mean(positions, axis=0) +[Lx/2,Ly/2,Lz/2]
    
    
    if invert:
        positions = positions[::-1]

    # Write new atom positions in the grofile
    new_lines = []
    n=0
    with open(grofile, 'r') as f:
        lines = f.readlines()

        for line in lines:
            words = line.split()
            try:
                float(words[-1])
            except:
                new_lines += [line]
                continue
            if len(words)<4:
                new_lines += [line]
                continue
            
            index = int(line[-30:-25])
            if (index < start_from_nth_atom+1) or (n>=len(positions)):
                new_lines += [line]
            else:
                line_ = line[:-25]+'%8.3f%8.3f%8.3f\n'%tuple(positions[n])
                new_lines += [line_]
                n=n+1

    data = ''.join(new_lines)
    with open('PA_box.gro', 'w') as f:    
        f.write(data)



def init_fiber_config_co(grofile, num_atomss, num_moleculess, Lx, Ly, Lz, start_from_nth_atoms, inverts, C_indicess,
    delta = 0.5, radial_offset = 0.5):
    """ Change the random positioning in grofile with fiber like initial config
    assumes all molecules are in contiguous in the grofile
    Note: start_from_nth_atom starts from index 0, whereas in grofile atom index starts from 1
    Assumes first residue is C
    C_indices: Carbon indices to use when aligning. Indices are after invertion if invert=True
    Vector C_indices[1]-C_indices[0] is used.
    If invert is true, the other end is pointed towards fiber inside

    num_atoms, mol_start_atom_indices, invert, C_indices are lists/arrays corresponding to grofiles' molecules
    which_mol are indices of molecules in grofiles corresponding to mol_start_atom_indices
    mol_indices are molecule number in the grofile

    delta = 0.5 # distance between layers in nm
    radial_offset = 0.5 # in nm

    NOTE: VERY ARDOUS CODE. NOT HAPPY
    """


    # Parameters for the initial configuration as a fiber 
    num_PA_layer = 9
    theta = np.pi*40/180 # angle between two PA in a layer
    theta_offset = np.pi*20/180
    

    atoms_positionss = []
    C_idss = []  # C_ids contain atom number of C in atom positions
    for k,num_atoms in enumerate(num_atomss):
        # Molecule's atoms' positions. Read from grofile
        atoms_positions = []
        invert = inverts[k]
        C_indices = C_indicess[k]
        start_from_nth_atom = start_from_nth_atoms[k]
        C_ids = []
        with open(grofile, 'r') as f:
            lines = f.readlines()
            start = False
            C_pos = []
            for line in lines[2:]:
                words = line.split()
                try:
                    float(words[-1])
                except:
                    continue
                if len(words)<4:
                    continue

                index = int(line[-30:-25])
                if index == start_from_nth_atom+1:
                    start = True
                if start:
                    if len(atoms_positions) == num_atoms:
                        break
                    
                    atoms_positions += [ [float(words[-3]), float(words[-2]), float(words[-1])] ]
                    if 'C' == words[1][0]:
                        C_pos += [ atoms_positions[-1] ]
                        C_ids += [len(atoms_positions)-1]
        if invert:
            atoms_positions = np.array(atoms_positions)#[::-1]
            C_ids = C_ids[::-1]
        else:
            atoms_positions = np.array(atoms_positions)
        

        # Translate with origin at outermost (C16) Carbon of alkyl tail
        atoms_positions -= atoms_positions[C_ids[C_indices[0]]]

        
        # Align atoms_positions in the x-y plane by aligning the vector
        # between first two consecutive C atoms of C16 towards x-axis
        # NOTE: Assumes the first atoms is C16
        v = np.array(atoms_positions[C_ids[C_indices[1]]]) - np.array(atoms_positions[C_ids[C_indices[0]]])
        q = quaternion.q_between_vectors(v, [1,0,0])
        for i,pos in enumerate(atoms_positions):
            atoms_positions[i] = quaternion.qv_mult(q, pos)
            
        atoms_positionss += [atoms_positions]
        C_idss += [C_ids] 
    
        
    C_positions = []  # of the outermost C  or innermost (of the core) C
    quaternions = [] # wrt when PAM chain is oriented in x-axis
    for i in range(sum(num_moleculess)):
        j = i  % num_PA_layer
        k = i // num_PA_layer
        
        q = quaternion.axisangle_to_q([0,0,1], theta * j + theta_offset * (1+(-1)**k)/2 )
        quaternions += [q]

        pos = np.array([0,0,k*delta]) + quaternion.qv_mult(q, [radial_offset,0,0])
        C_positions += [ pos ]
    

    # randomly arrange molecules
    args = np.arange(sum(num_moleculess))
    np.random.shuffle(args)
    bins = []
    for arg in args:
        s = 0
        for i,num_atoms in enumerate(num_atomss):
            if arg < s + num_moleculess[i]:
                bins+=[i]
                break
            else:
                s += num_moleculess[i]
    

    # Calculate positions of all atoms
    t = np.array(num_atoms)*np.array(num_moleculess)
    positions = []
    for i,arg in enumerate(args):
        atoms_positions = atoms_positionss[bins[i]]
        invert = inverts[bins[i]]
        positions_ = []
        for j,pos in enumerate(atoms_positions):
            p = quaternion.qv_mult(quaternions[i], pos)
            p += C_positions[i]
            positions_ += [ p ]
        # if invert:
        #     positions += positions_#[::-1]
        # else:
        #     positions += positions_
        positions += [ positions_ ]


    positions_ = [[]]*len(positions)
    for i,pos in enumerate(positions):
        positions_[args[i]] = pos
    positions = positions_
    
    positions_=[]
    for pos in positions:
        positions_ += pos
    positions = np.array(positions_)
    

    positions = np.array(positions) - np.mean(positions, axis=0) +[Lx/2,Ly/2,Lz/2]
    

    # Write new atom positions in the grofile
    new_lines = []
    n=0
    start_from_nth_atom = 0 
    with open(grofile, 'r') as f:
        lines = f.readlines()

        for line in lines:
            words = line.split()
            try:
                float(words[-1])
            except:
                new_lines += [line]
                continue
            if len(words)<4:
                new_lines += [line]
                continue
            
            index = int(line[-30:-25])
            if (index < start_from_nth_atom+1) or (n>=len(positions)):
                new_lines += [line]
            else:
                line_ = line[:-25]+'%8.3f%8.3f%8.3f\n'%tuple(positions[n])
                new_lines += [line_]
                n=n+1

    data = ''.join(new_lines)
    with open(grofile, 'w') as f:    
        f.write(data)

        
def init_lamella_config(grofile, num_atoms, num_molecules, Lx, Ly, Lz, start_from_nth_atom=0, invert=False, C_indices=[0,2]):
    """ Change the random positioning in grofile with lamella membrane like initial config
    assumes all molecules are in contiguous in the grofile.
    All molecules are arranged parallel to each other
    Note: start_from_nth_atom starts from index 0, whereas in grofile atom index starts from 1
    Assumes first residue is C
    C_indices: Carbon indices to use when aligning. 
    Vector C_indices[1]-C_indices[0] is used.
    If invert is true, the other end is pointed towards fiber inside
    """


    # Parameters for the initial configuration as a fiber 
    num_PA_layer = int(np.sqrt(num_molecules))
    # theta = np.pi*40/180 # angle between two PA in a layer
    # theta_offset = np.pi*20/180
    delta = 0.4 # distance between molecules in nm
    # radial_offset = 0.5 # in nm


    # Molecule's atoms' positions. Read from grofile
    atoms_positions = []
    with open(grofile, 'r') as f:
        lines = f.readlines()
        start = False
        C_pos = []
        for line in lines[2:]:
            words = line.split()
            try:
                float(words[-1])
            except:
                continue
            if len(words)<4:
                continue

            index = int(line[-30:-25])
            if index == start_from_nth_atom+1:
                start = True
            if start:
                if len(atoms_positions) == num_atoms:
                    break
                atoms_positions += [ [float(words[-3]), float(words[-2]), float(words[-1])] ]
                if 'C' == words[1][0]:
                    C_pos += [ atoms_positions[-1] ]
    
    if invert:
        atoms_positions = np.array(atoms_positions)[::-1]
    else:
        atoms_positions = np.array(atoms_positions)

    # Translate with origin at outermost (C16) Carbon of alkyl tail
    atoms_positions -= atoms_positions[C_indices[0]]

    # Align atoms_positions in the x-y plane by aligning the vector 
    # between first two consecutive C atoms of C16 towards x-axis
    # NOTE: Assumes the first atoms is C16
    v = np.array(C_pos[C_indices[1]]) - np.array(C_pos[C_indices[0]])
    q = quaternion.q_between_vectors(v, [0,0,1])
    for i,pos in enumerate(atoms_positions):
        atoms_positions[i] = quaternion.qv_mult(q, pos)
    
        
    C_positions = []  # of the outermost C
    quaternions = [] # wrt when PAM chain is oriented in x-axis
    for i in range(num_molecules):
        j = i  % num_PA_layer
        k = i // num_PA_layer
        
        q = [1,0,0,0]
        quaternions += [q]

        pos = np.array([j*delta, k*delta, 0]) #+ quaternion.qv_mult(q, [radial_offset,0,0])
        C_positions += [ pos ]
        
    
    # Calculate positions of all atoms
    positions = []
    for i in range(num_molecules):
        for j,pos in enumerate(atoms_positions):
            p = quaternion.qv_mult(quaternions[i], pos)
            p += C_positions[i]
            positions += [ p ]
    positions = np.array(positions) - np.mean(positions, axis=0) +[Lx/2,Ly/2,Lz/2]
    
    
    if invert:
        positions = positions[::-1]

    # Write new atom positions in the grofile
    new_lines = []
    n=0
    with open(grofile, 'r') as f:
        lines = f.readlines()

        for line in lines:
            words = line.split()
            try:
                float(words[-1])
            except:
                new_lines += [line]
                continue
            if len(words)<4:
                new_lines += [line]
                continue
            
            index = int(line[-30:-25])
            if (index < start_from_nth_atom+1) or (n>=len(positions)):
                new_lines += [line]
            else:
                line_ = line[:-25]+'%8.3f%8.3f%8.3f\n'%tuple(positions[n])
                new_lines += [line_]
                n=n+1

    data = ''.join(new_lines)
    with open('PA_box.gro', 'w') as f:    
        f.write(data)
        

def init_bilayer_config(grofile, num_atoms, num_molecules, Lx, Ly, Lz, start_from_nth_atom=0, invert=False, C_indices=[0,2]):
    """ Change the random positioning in grofile with lamella membrane like initial config
    assumes all molecules are in contiguous in the grofile.
    All molecules are arranged parallel to each other
    Note: start_from_nth_atom starts from index 0, whereas in grofile atom index starts from 1
    Assumes first residue is C
    C_indices: Carbon indices to use when aligning. 
    Vector C_indices[1]-C_indices[0] is used.
    If invert is true, the other end is pointed towards fiber inside
    """
    use_one_less = False
    if num_molecules%2!=0:
        use_one_less = True
        num_molecules = int(num_molecules/2)+1
    else:
        num_molecules = int(num_molecules/2)
    

    # Parameters for the initial configuration as a fiber 
    num_PA_layer = int(np.sqrt(num_molecules))
    # theta = np.pi*40/180 # angle between two PA in a layer
    # theta_offset = np.pi*20/180
    delta = 0.4 # distance between molecules in nm
    # radial_offset = 0.5 # in nm


    # Molecule's atoms' positions. Read from grofile
    atoms_positions = []
    with open(grofile, 'r') as f:
        lines = f.readlines()
        start = False
        C_pos = []
        for line in lines[2:]:
            words = line.split()
            try:
                float(words[-1])
            except:
                continue
            if len(words)<4:
                continue

            index = int(line[-30:-25])
            if index == start_from_nth_atom+1:
                start = True
            if start:
                if len(atoms_positions) == num_atoms:
                    break
                atoms_positions += [ [float(words[-3]), float(words[-2]), float(words[-1])] ]
                if 'C' == words[1][0]:
                    C_pos += [ atoms_positions[-1] ]
    
    if invert:
        atoms_positions = np.array(atoms_positions)[::-1]
    else:
        atoms_positions = np.array(atoms_positions)

    # Translate with origin at outermost (C16) Carbon of alkyl tail
    atoms_positions -= atoms_positions[C_indices[0]]

    # Align atoms_positions in the x-y plane by aligning the vector 
    # between first two consecutive C atoms of C16 towards x-axis
    # NOTE: Assumes the first atoms is C16
    v = np.array(C_pos[C_indices[1]]) - np.array(C_pos[C_indices[0]])
    q = quaternion.q_between_vectors(v, [0,0,1])
    for i,pos in enumerate(atoms_positions):
        atoms_positions[i] = quaternion.qv_mult(q, pos)
    
        
    C_positions = []  # of the outermost C
    quaternions = [] # wrt when PAM chain is oriented in x-axis
    for i in range(num_molecules):
        j = i  % num_PA_layer
        k = i // num_PA_layer
        
        q = [1,0,0,0]
        quaternions += [q]

        pos = np.array([j*delta, k*delta, 0]) #+ quaternion.qv_mult(q, [radial_offset,0,0])
        C_positions += [ pos ]


    # Calculate positions of all atoms
    positions = []
    for i in range(num_molecules):
        for j,pos in enumerate(atoms_positions):
            p = quaternion.qv_mult(quaternions[i], pos)
            p += C_positions[i]
            positions += [ p ]
    
    # add bottom layer of bilayer
    q = quaternion.axisangle_to_q([0,1,0], np.pi )
    positions_=[]
    for i in range(num_molecules):
        for j in range(num_atoms):
            pos = np.copy(positions[i*num_atoms+j])
            positions_ += [ np.array([0,0,-0.1]) + C_positions[i] + quaternion.qv_mult(q, pos-C_positions[i]) ]
    
    if use_one_less:
        positions += positions_[:-1*num_atoms]
    else:
        positions += positions_

    positions = np.array(positions) - np.mean(positions, axis=0) +[Lx/2,Ly/2,Lz/2]
    

    if invert:
        positions = positions[::-1]

    # Write new atom positions in the grofile
    new_lines = []
    n=0
    with open(grofile, 'r') as f:
        lines = f.readlines()

        for line in lines:
            words = line.split()
            try:
                float(words[-1])
            except:
                new_lines += [line]
                continue
            if len(words)<4:
                new_lines += [line]
                continue
            
            index = int(line[-30:-25])
            if (index < start_from_nth_atom+1) or (n>=len(positions)):
                new_lines += [line]
            else:
                line_ = line[:-25]+'%8.3f%8.3f%8.3f\n'%tuple(positions[n])
                new_lines += [line_]
                n=n+1

    data = ''.join(new_lines)
    with open('PA_box.gro', 'w') as f:    
        f.write(data)
    
    
def get_num_atoms_fromGRO(filename):
    # Returns the number in the 2nd line of grofile

    num_atoms = 0
    
    with open(filename, 'r') as f:
        lines = f.readlines()
        num_atoms = int(lines[1].split()[0])

    return num_atoms


def get_num_molecules_fromtop(topfile, molname):
    # Get the number of Protein from the [ molecules ] section

    with open(topfile, 'r') as f:
        lines = f.readlines()
        start = False
        for i,line in enumerate(lines):
            if '#mols' in line:
                words = lines[i+1].split()
                if molname == words[0]:
                    num_molecules = int(words[1])
                    break
        return num_molecules


################################################################################
