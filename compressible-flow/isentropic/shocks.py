#!/usr/bin/env python

# Fully calculates post normal shock state given Mach and gamma
# Ratios can be calculated without conditions
# Post shock conditions require pre shock conditions

# Does not check for errors
# Does not solve for pre-shock mach number numerically, requires pressure ratio
# To Do: 
# should be able to only give one of the conditions and gamma and solve for everything else

class NormalShock:

    """
    Normal shock class needs Mach, gamma input to calculate property ratios
    Given pre-shock conditions, can get post-shock conditions
    """

    def __init__(self,M1=None,gamma=1.4,p1_static=None,rho1_static=None,T1_static=None,p1_total=None,T1_total=None,rho1_total=None):
        self.gamma = gamma
        self.M1 = M1
        self.M2 = get_mach_normal_shock(M1=self.M1,gamma=self.gamma)
        self.p1_static = p1_static
        self.p2_static = get_static_pressure_normal_shock(M1=self.M1,gamma=self.gamma,p1_static=self.p1_static)
        self.p_ratio_static = get_static_pressure_ratio_normal_shock(M1=self.M1,gamma=self.gamma)
        self.rho1_static = rho1_static
        self.rho2_static = get_static_density_normal_shock(M1=self.M1,gamma=self.gamma,rho1_static=self.rho1_static)
        self.rho_ratio_static = get_static_density_ratio_normal_shock(M1=self.M1,gamma=self.gamma)
        self.T1_static = T1_static
        self.T2_static = get_static_temperature_normal_shock(M1=self.M1,gamma=self.gamma,T1_static=self.T1_static)
        self.T_ratio_static = get_static_temperature_ratio_normal_shock(M1=self.M1,gamma=self.gamma)
        self.p1_total = p1_total
        self.p2_total = get_total_pressure_normal_shock(M1=self.M1,gamma=self.gamma,p1_total=self.p1_total)
        self.p_ratio_total = get_total_pressure_ratio_normal_shock(M1=self.M1,gamma=self.gamma)
        self.T1_total = T1_total
        self.T2_total = get_total_temperature_normal_shock(T1_total=self.T1_total)
        self.T_ratio_total = get_total_temperature_ratio_normal_shock()
        self.rho1_total = rho1_total
        self.rho2_total = get_total_density_normal_shock(p_ratio_total=self.p_ratio_total)
        self.rho_ratio_total = get_total_density_ratio_normal_shock(p_ratio_total=self.p_ratio_total)
        self.sonic_area_ratio = get_sonic_area_ratio_normal_shock(p_ratio_total=self.p_ratio_total)
        self.p2_total_p1_static = get_p2_total_over_p1_static(M2=self.M2,gamma=self.gamma,p_ratio_static=self.p_ratio_static)

def get_mach_normal_shock(M1=None,gamma=1.4):
    M2 = ((1 + (gamma-1)/2 * M1**2) / (gamma * M1**2 - (gamma-1)/2))**0.5
    print(f'Post normal-shock Mach number = {M2}') 
    return M2

def get_upstream_mach_normal_shock(M2=None,gamma=1.4,p_ratio_static=None):
    M1 = ((p_ratio_static + (gamma-1)/(gamma+1))*(gamma+1)/(2*gamma))**0.5
    print(f'Pre shock mach number = {M1}')
    return M1

def get_static_pressure_normal_shock(M1=None,gamma=1.4,p1_static=None):

    if p1_static:

        p2_static = (1 + 2*gamma/(gamma+1)*(M1**2-1))*p1_static
        print(f'Post shock static pressure = {p2_static}')
        return p2_static

def get_static_pressure_ratio_normal_shock(M1=None,gamma=1.4):
    p_ratio_static = (1 + 2*gamma/(gamma+1)*(M1**2-1))
    print(f'Post shock static pressure ratio = {p_ratio_static}')
    return p_ratio_static

def get_static_density_normal_shock(M1=None,gamma=1.4,rho1_static=None):

    if rho1_static:

        rho2_static = (gamma+1)*M1**2/(2+(gamma-1)*M1**2)*rho1_static
        print(f'Post shock static density = {rho2_static}')
        return rho2_static

def get_static_density_ratio_normal_shock(M1=None,gamma=1.4):
    rho_ratio_static = (gamma+1)*M1**2/(2+(gamma-1)*M1**2)
    print(f'Post shock static density ratio = {rho_ratio_static}')
    return rho_ratio_static

def get_static_temperature_normal_shock(M1=None,gamma=1.4,T1_static=None):

    if T1_static:

        T2_static = (1 + 2*gamma/(gamma+1) * (M1**2-1)) *  ((2 + (gamma-1) * M1**2)/((gamma+1)*M1**2))*T1_static
        print(f'Post shock static temperature = {T2_static}')
        return T2_static

def get_static_temperature_ratio_normal_shock(M1=None,gamma=1.4):
    T_ratio_static = (1 + 2*gamma/(gamma+1) * (M1**2-1)) *  ((2 + (gamma-1) * M1**2)/((gamma+1)*M1**2))
    print(f'Post shock static temperature ratio = {T_ratio_static}')
    return T_ratio_static

def get_total_pressure_normal_shock(M1=None,gamma=1.4,p1_total=None):

    if p1_total:

        p2_total = (((gamma+1)/2*M1**2)/(1 + (gamma-1)/2 * M1**2))**(gamma/(gamma-1)) * (2*gamma/(gamma+1)*M1**2-(gamma-1)/(gamma+1))**(1/(1-gamma))*p1_total
        print(f'Post shock total pressure ratio = {p2_total}')
        return p2_total    

def get_total_pressure_ratio_normal_shock(M1=None,gamma=1.4):
    p_ratio_total = (((gamma+1)/2*M1**2)/(1 + (gamma-1)/2 * M1**2))**(gamma/(gamma-1)) * (2*gamma/(gamma+1)*M1**2-(gamma-1)/(gamma+1))**(1/(1-gamma))
    print(f'Post shock total pressure ratio = {p_ratio_total}')
    return p_ratio_total

def get_total_temperature_normal_shock(T1_total=None):

    if T1_total:

        T2_total = T1_total
        print(f'Post shock total temperature = {T2_total}')
        return T2_total

def get_total_temperature_ratio_normal_shock():
    T_ratio_total = 1
    print(f'Post shock total temperature ratio = {T_ratio_total}')
    return T_ratio_total

def get_total_density_normal_shock(p_ratio_total=None,rho1_total=None):

    if rho1_total:

        rho2_total = p_ratio_total*rho1_total
        print(f'Post shock total density = {rho2_total}')
        return rho2_total

def get_total_density_ratio_normal_shock(p_ratio_total=None):
    rho_ratio_total = p_ratio_total
    print(f'Post shock total density ratio = {rho_ratio_total}')
    return rho_ratio_total

def get_sonic_area_ratio_normal_shock(p_ratio_total=None):
    sonic_area_ratio = 1/p_ratio_total
    print(f'Post shock sonic area ratio = {sonic_area_ratio}')
    return sonic_area_ratio

def get_p2_total_over_p1_static(M2=None,gamma=1.4,p_ratio_static=None):
    p2_total_p1_static = (1 + (gamma-1)/2*M2**2)**(gamma/(gamma-1))*p_ratio_static
    print(f'P2t/P1 = {p2_total_p1_static}')
    return p2_total_p1_static


if __name__=='__main__':
    #mach = get_mach_normal_shock(M1=1.5,gamma=1.4)
    #p2 = get_static_pressure_normal_shock(M1=1.5,p1_static=101322)
    #p2r = get_static_pressure_ratio_normal_shock(M1=1.5)
    #rho2 = get_static_density_normal_shock(M1=1.5,rho1_static=1.22)
    #rho2r = get_static_density_ratio_normal_shock(M1=1.5)
    #T2 = get_static_temperature_normal_shock(M1=1.5,T1_static=300)
    #T2r = get_static_temperature_ratio_normal_shock(M1=1.5)
    #ptr = get_total_pressure_ratio_normal_shock(M1=1.5)
    #p2t_p1 = get_p2_total_over_p1_static(M2=mach,p_ratio_static=p2r)
    #M1 = get_upstream_mach_normal_shock(M2=0.8,p_ratio_static=1.72)

    foo = NormalShock(M1=1.5)