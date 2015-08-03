from enum import Enum;

class GEOMETRY:
	SPHERICAL   = 'SPHERICAL';
	CYLINDRICAL = 'CYLINDRICAL';
	CARTESIAN = 'CARTESIAN';

class DENSITY:
	NORMAL = 'DENSITY_NORMAL';
	BOUSSINESQ = 'DENSITY_BOUSSINESQ';
	
class OPERATOR_COMPLEXITY:
	BARO_ALL = 'BARO_ALL';
	BAROTROPIC = 'BAROTROPIC';
	BAROCLINIC = 'BAROCLINIC';