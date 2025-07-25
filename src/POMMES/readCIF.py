# PARAMETERS
disordered_file='../../data/EuTaO3disordered.cif'
#CIF file of ordered structure
ordered_file='../../data/EuTaO3ordered.cif'
import numpy as np
def filename(path):
    """Extract filename from filepath"""
    try:
        return path.split("/")[-1]
    except:
        return ''
def warn_user(alert,prompt):
    """Print warning and ask user whether or not to continue despite that"""
    print("Warning:",alert)
    if (input(prompt+' y/[n]'+'\t')[0].lower()!='y'):
            raise Exception('Exit.')
def line_index_of(lines,string,start=0,file=''):
    """
    Return index of first line in an list of lines that contains `string`
    
    Parameters
    ----------
    lines : list[str]
        Lines to search through.
    string : str
        String to look for.
    start : int, optional
        Line index to start searching at, default 0.
    file : str, optional
        Descriptive name.    

    Raises
    ------
    IndexError
        If no line contains `string`.

    Returns
    -------
    i : int
        Index in `lines` of the first line after `lines[start]` containing `string`.
    """
    if (len(string)<1):
        raise IndexError('Cannot search for empty string')
    length=len(lines)
    i=start
    if (lines[i].find(string)>-1):
        return i
    i=(i+1)%length
    while (i!=start):
        if (lines[i].find(string)>-1):
            return i
        else:
            i=(i+1)%length
    raise IndexError('String '+string+' not found in file '+file)
def validate_spacegroup(filename,lines):
    """
    Find the space group in a file and check that it is P 1.

    Parameters
    ----------
    filename : str
        Name of file whose space group is being checked.
    lines : list[str]
        List of lines in the file.

    Raises
    ------
    Exception
        Quit if the user does not want to continue.

    Returns
    -------
    spacegroupline : int
        The line index at which the space group is found.

    """
    spacegroupline=0
    try:
        spacegroupline=line_index_of(lines,'space_group_name',file=filename)
        if (lines[spacegroupline].find('P 1')<0):
            warn_user('''It seems the file '''+filename+''' is not asymmetrized (space group P 1)
    This code requires P 1 in order to have all atomic positions explicitly shown in the CIF file.
    Note: certain crystal visualization programs let you remove symmetry without losing atoms in order to achieve this.''','Continue?')
    except IndexError:
        warn_user('No space group seems to be specified in the file '+filename+'.', 'Continue, assuming P 1?')
    return spacegroupline
def get_lattice_constants(path,is_path=True):
    """
    Get the lattice constants a,b,c of a CIF file.

    Parameters
    ----------
    path : str
        Path to a CIF file to read.
    is_path : boolean
        If True, `path` is a filepath, else it is the file's contents given by readlines().
        
    Returns
    -------
    list[float]
        Latice constants [a,b,c].

    """
    lines=None
    if is_path:
        cif=open(path,'r')
        lines=cif.readlines()
        cif.close()
    else:
        lines=path
    p=path if is_path else ''
    a=line_index_of(lines, "_cell_length_a",file=p)
    b=line_index_of(lines, "_cell_length_b",file=p,start=a)
    c=line_index_of(lines, "_cell_length_c",file=p,start=b)
    return list(map(lambda i: float(lines[i].split()[1].replace("(","").replace(")","")),[a,b,c]))

def make_cif_matrices(disordered_file,ordered_file,is_path=True):
    """
    Return matrices for the sites and atomic positions/occupancies in the
    disordered and ordered structures for interpolation purposes.
    
    If a site/atom combination appears in both files, both its position and
    occupancy will interpolate between the two.
    
    If a site/atom combination appears only in one file, it is interpolated to
    a vacancy (0 occupancy) that is added to the other file's matrix in the
    same position.

    Parameters
    ----------
    disordered_file : str
        Path to a CIF file for the disordered structure.
    ordered_file : str
        Path to a CIF file for the ordered structure.
    is_path : boolean
        If True, the above files are paths, else they are the files' contents from readlines().

    Returns
    -------
    A tuple of the following 3 values:
        
    sites : np.ndarray
        Matrix where column 0 is the list of atom site labels in the
        order they appear in the atom matrices and column 1 is the atom/ion
        types in that same order.
    disordered : np.ndarray
        Matrix where columns 0, 1, and 2 are the x, y, and z fractional
        coordinates, and column 3 is the occupancy, of each unique site/atom
        combination in the disordered CIF in the order given by `sites`.
    ordered : np.ndarray
        Same as `disordered` but for the ordered CIF.

    """
    linesD=None
    linesO=None
    if is_path:
        disorderedCIF=open(disordered_file,'r') if is_path else disordered_file
        orderedCIF=open(ordered_file,'r') if is_path else ordered_file
        linesD=disorderedCIF.readlines()
        linesO=orderedCIF.readlines()
        disorderedCIF.close()
        orderedCIF.close()
    else:
        linesD=disordered_file
        linesO=ordered_file
    sites=[]
    def read_atoms(lines,filename='',start=0):
        pos=line_index_of(lines,'atom_site_label',start=start,file=filename)
        while (lines[pos].find('loop_')>=0):
            pos-=1
        col_index=-1
        label_col=0
        atom_col=1
        x_col=4
        y_col=5
        z_col=6
        occ_col=8
        #Find which columns of the atom "table" contain the needed information
        #since they are not always in the same order
        while lines[pos].strip()[0]=='_':
            col_index+=1
            if (lines[pos].find("atom_site_label")>0):
                label_col=col_index
            elif (lines[pos].find("atom_site_type_symbol")>0):
                atom_col=col_index
            elif (lines[pos].find("atom_site_occupancy")>0):
                occ_col=col_index
            elif (lines[pos].find("atom_site_fract_x")>0):
                x_col=col_index
            elif (lines[pos].find("atom_site_fract_y")>0):
                y_col=col_index
            elif (lines[pos].find("atom_site_fract_z")>0):
                z_col=col_index
            pos+=1
        atoms=[]
        for l in range(pos,len(lines)):
            if (lines[l].strip()=='' or lines[l][0]=='#'):
                break
            line=lines[l].split()
            label=line[label_col]
            atom=line[atom_col]
            coordinates=tuple(map(lambda v: float(v.replace("(","").replace(")","")),[line[x_col],line[y_col],line[z_col]]))
            occupancy=float(line[occ_col].replace("(","").replace(")",""))
            # print("Site:",label,"Atom:",atom,"coordinates:",coordinates,"occupancy:",occupancy)
            atoms.append([label,atom,coordinates[0],coordinates[1],coordinates[2],occupancy])
            if (sites.count([label,atom])<1):
                sites.append([label,atom])
        return atoms
    atomsD=read_atoms(linesD,filename(disordered_file),validate_spacegroup(filename(disordered_file),linesD))
    atomsO=read_atoms(linesO,filename(ordered_file),validate_spacegroup(filename(ordered_file),linesO))
    ordered=[]
    disordered=[]
    for s in sites:
        matchD=False
        for d in atomsD:
            if d[0:2]==s:
                matchD=True
                rowD=d[2:]
                break
        if not matchD:
            print('not matchD')
            rowD=[0,0,0,0]
        matchO=False
        for o in atomsO:
            if o[0:2]==s:
                matchO=True
                rowO=o[2:]
                break
        if matchO:
            if not matchD:
                #Site only in ordered structure, assume vacancy in disordered
                rowD[0:3]=rowO[0:3]
        else:
            if matchD:
                #Site only in disordered structure, assume vacancy in ordered
                rowO=rowD[0:3]+[0]
            else:
                #Site somehow absent in both
                print('If you see this, something is wrong with your matrices.')
                rowO=[0,0,0,0]
        disordered.append(rowD)
        ordered.append(rowO)
    return np.array(sites), np.array(disordered), np.array(ordered)

if __name__ == "__main__":
    sites, disordered, ordered = make_cif_matrices(disordered_file,ordered_file)
    print('ATOMS\t\t\t\tDISORDERED\t\t\t\tORDERED')
    for s in range(len(sites)):
        print(sites[s],end='\t\t')
        print(disordered[s],end='\t\t')
        print(ordered[s],end='\n')
    print('\nPreload code:')
    print("import numpy as np\nsites=np.array(",sites.tolist(),")\ndisordered=np.array(",disordered.tolist(),")\nordered=np.array(",ordered.tolist(),")",sep="")
    
    '''disordered and ordered are matrices of the disordered and ordered atoms'
    coordinates and occupancies. Atoms only present in one are interpolated to
    vacancies in the same positions. Atoms present in both have their positions
    and occupancies interpolated.'''