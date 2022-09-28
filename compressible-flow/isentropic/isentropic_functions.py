#!/usr/bin/env python

from shocks import NormalShock

# Currently will run in a loop until all states fully defined
# Does not check for compatibility between states
# Does not catch errors
# Need to add option to input Mach and gamma to get total/static ratios without looping forever

# To Do:
# Add sonic ratios

class CompressibleFlow:

    def __init__(self,M=None,u=None,a=None,p=None,rho=None,T=None,p_t=None,rho_t=None,T_t=None,p_t_ratio=None,rho_t_ratio=None,T_t_ratio=None,gamma=1.4,gas='Air',R=287):
        self.M = M
        self.u = u
        self.a = a
        self.p = p
        self.rho = rho
        self.T = T
        self.p_t = p_t
        self.rho_t = rho_t
        self.T_t = T_t
        self.p_t_ratio = p_t_ratio
        self.rho_t_ratio = rho_t_ratio
        self.T_t_ratio = T_t_ratio
        self.gamma = gamma
        self.gas = gas
        self.R = R

    def isentropic_state(self):

        mode = self.__isentropic_input_check()
        print(f'Calculation mode: {mode}')

        if mode == 0:
            self.__isentropic_state_method_0()
                
        elif mode == 1:
            self.__isentropic_state_method_1()

        elif mode == 2:
            self.__isentropic_state_method_2()

    def __isentropic_state_method_0(self):
        self.p_t_ratio = get_pressure_ratio(M=self.M,gamma=self.gamma,R=self.R)
        self.rho_t_ratio = get_density_ratio(M=self.M,gamma=self.gamma,R=self.R)
        self.T_t_ratio = get_temperature_ratio(M=self.M,gamma=self.gamma,R=self.R)
        
    def __isentropic_state_method_1(self):

        while None in self.__dict__.values():
            if not self.p:
                self.p = get_static_pressure(gamma=self.gamma,R=self.R,M=self.M,u=self.u,rho=self.rho,T=self.T,p_t=self.p_t,rho_t=self.rho_t,T_t=self.T_t,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.rho:
                self.rho = get_static_density(gamma=self.gamma,R=self.R,M=self.M,u=self.u,p=self.p,T=self.T,p_t=self.p_t,rho_t=self.rho_t,T_t=self.T_t,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.T:
                self.T = get_static_temperature(gamma=self.gamma,R=self.R,M=self.M,u=self.u,p=self.p,rho=self.rho,p_t=self.p_t,rho_t=self.rho_t,T_t=self.T_t,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.a:
                self.a = get_sonic_velocity(gamma=1.4,R=287,u=self.u,p=self.p,rho=self.rho,T=self.T,p_t=self.p_t,rho_t=self.rho_t,T_t=self.T_t)
            if not self.M:
                self.M = get_mach_number(gamma=1.4,R=287,a=self.a,u=self.u,p=self.p,rho=self.rho,T=self.T,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.u:
                self.u = get_fluid_velocity(gamma=1.4,R=287,M=self.M,a=self.a,p=self.p,rho=self.rho,T=self.T)
            if not self.p_t:
                self.p_t = get_total_pressure(gamma=self.gamma,R=self.R,M=self.M,u=self.u,p=self.p,rho=self.rho,T=self.T,rho_t=self.rho_t,T_t=self.T_t,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.rho_t:
                self.rho_t = get_total_density(gamma=self.gamma,R=self.R,M=self.M,u=self.u,p=self.p,rho=self.rho,T=self.T,p_t=self.p_t,T_t=self.T_t,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.T_t:
                self.T_t = get_total_temperature(gamma=self.gamma,R=self.R,M=self.M,u=self.u,p=self.p,rho=self.rho,T=self.T,p_t=self.p_t,rho_t=self.rho_t,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.p_t_ratio:
                self.p_t_ratio = get_pressure_ratio(gamma=self.gamma,R=self.R,M=self.M,u=self.u,p=self.p,rho=self.rho,T=self.T,p_t=self.p_t,rho_t=self.rho_t,T_t=self.T_t,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.rho_t_ratio:
                self.rho_t_ratio = get_density_ratio(gamma=self.gamma,R=self.R,M=self.M,u=self.u,p=self.p,rho=self.rho,T=self.T,p_t=self.p_t,rho_t=self.rho_t,T_t=self.T_t,p_t_ratio=self.p_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.T_t_ratio:
                self.T_t_ratio = get_temperature_ratio(gamma=self.gamma,R=self.R,M=self.M,u=self.u,p=self.p,rho=self.rho,T=self.T,p_t=self.p_t,rho_t=self.rho_t,T_t=self.T_t,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio)

    def __isentropic_state_method_2(self):
        while None in [self.M,self.p_t_ratio,self.rho_t_ratio,self.T_t_ratio]:
            if not self.M:
                self.M = get_mach_number(gamma=self.gamma,p_t_ratio=self.p_t_ratio,rho_t_ratio=self.rho_t_ratio,T_t_ratio=self.T_t_ratio)
            if not self.p_t_ratio:
                self.p_t_ratio = get_pressure_ratio(M=self.M,gamma=self.gamma)
            if not self.rho_t_ratio:
                self.rho_t_ratio = get_density_ratio(M=self.M,gamma=self.gamma)
            if not self.T_t_ratio:
                self.T_t_ratio = get_temperature_ratio(M=self.M,gamma=self.gamma)
            pass

    def __isentropic_input_check(self):

        mode = None

        if self.M:

            if not any([self.p,self.rho,self.T]):
                mode = 0
                print('Ideal gas not fully defined -- calculating thermodynamic ratios only.')
                return mode
            
            elif sum(1 for i in [self.p,self.rho,self.T] if i is not None)>1:
                mode = 1
                print('Ideal gas fully defined -- calculating all values.')
                return mode

        elif sum(1 for i in [self.p_t_ratio,self.rho_t_ratio,self.T_t_ratio] if i is not None)>0:
            mode = 2
            print('No Mach number defined -- using given total ratio to calculate other thermodynamic ratios.')
            return mode

        elif not self.M:
            raise ValueError('No Mach number defined, ideal gas not fully defined, no total ratios defined -- exiting now.')
        
    def report_outputs(self):
        for key,value in zip(self.__dict__.keys(),self.__dict__.values()):
            print(str(key) + ': ' + str(value))

    def shock(self):
        shocked = NormalShock(M1=self.M,gamma=self.gamma,p1_static=self.p,rho1_static=self.rho_static,T1_static=self.T_static,p1_total=self.p_total,T1_total=self.T_total,rho1_total=self.rho_total)
        post_shock = CompressibleFlow(M=shocked.M2,p=shocked.p2_static,rho_static=shocked.rho2_static,T_static=shocked.T2_static)
        return shocked, post_shock


def get_sonic_velocity(gamma=1.4,R=287,u=None,p=None,rho=None,T=None,p_t=None,rho_t=None,T_t=None):

    a = None
    
    while not a:

        if T:
            a = (gamma*R*T)**.5

        elif p and rho:
            a = (gamma*p/rho)**0.5

        elif u == 0 and T_t:
            a = (gamma*R*T_t)**.5
                
        elif u == 0 and p_t and rho_t:
            a = (gamma*p_t/rho_t)**0.5

    print(f'Speed of sound = {a}')
    return a

def get_fluid_velocity(gamma=1.4,R=287,M=None,a=None,p=None,rho=None,T=None):

    u = None

    while not u:

        if M and a:
            u = M*a

        elif M and not a and T:
            u = M*(gamma*R*T)**0.5

        elif M and not T and p and rho:
            u = M*(gamma*p/rho)**0.5

        elif M == 0:
            u = 0

    print(f'Fluid velocity = {u}')
    return u

def get_mach_number(gamma=1.4,R=287,a=None,u=None,p=None,rho=None,T=None,p_t_ratio=None,rho_t_ratio=None,T_t_ratio=None):
    """
    Calculates Mach number
    """
    M = None

    while not M:

        if u and a:
            M = u/a

        elif u and not a and T:
            M = u/(gamma*R*T)**0.5

        elif u and not a and not T and p and rho:
            M = u/(gamma*p/rho)**0.5

        elif p_t_ratio:
            M = ((p_t_ratio**((gamma-1)/gamma)-1)*2/(gamma-1))**0.5
            pass

        elif rho_t_ratio:
            M = ((rho_t_ratio**(gamma-1)-1)*2/(gamma-1))**0.5
            pass 

        elif T_t_ratio:
            M = ((T_t_ratio-1)*2/(gamma-1))**0.5
            pass

        elif u == 0:
            M = 0

    print(f'Mach number= {M}')
    return M

def get_static_pressure(gamma=1.4,R=287,M=None,u=None,rho=None,T=None,p_t=None,rho_t=None,T_t=None,p_t_ratio=None,rho_t_ratio=None,T_t_ratio=None):

    p = None

    while not p:

        if rho and T:
            p = rho*R*T
        
        elif p_t and M:
            p = p_t/(1 + (gamma-1) / 2 * M**2)**(gamma/(gamma-1))

        elif p_t_ratio and p_t:
            p = p_t/p_t_ratio

        elif p_t and T_t_ratio:
            p = p_t/(T_t_ratio**(gamma/(gamma-1)))

        elif p_t and T and T_t:
            p = p_t/((T_t/T)**(gamma/(gamma-1)))

        elif p_t and rho_t_ratio:
            p = p_t/rho_t_ratio**gamma

        elif p_t and rho and rho_t:
            p = p_t/(rho_t/rho)**gamma

        elif p_t and (u == 0):
            p = p_t

    print(f'Static pressure is {p}')
    return p

def get_total_pressure(gamma=1.4,R=287,M=None,u=None,p=None,rho=None,T=None,rho_t=None,T_t=None,p_t_ratio=None,rho_t_ratio=None,T_t_ratio=None):

    p_t = None

    while not p_t:

        if rho_t and T_t:
            p_t = rho_t*R*T_t

        elif p and M:
            p_t = (1 + (gamma-1) / 2 * M**2)**(gamma/(gamma-1))*p

        elif p_t_ratio and p:
            p_t = p_t_ratio*p

        elif p and T_t_ratio:
            p_t = p*(T_t_ratio**(gamma/(gamma-1)))

        elif p and T and T_t:
            p_t = p*((T_t/T)**(gamma/(gamma-1)))

        elif p and rho_t_ratio:
            p_t = p*rho_t_ratio**gamma

        elif p and rho and rho_t:
            p_t = p*(rho_t/rho)**gamma     

        elif p and (u == 0):
            p_t = p
    
    print(f'Total pressure = {p_t}')
    return p_t

def get_pressure_ratio(gamma=1.4,R=287,M=None,u=None,p=None,rho=None,T=None,p_t=None,rho_t=None,T_t=None,rho_t_ratio=None,T_t_ratio=None):

    p_t_ratio = None

    while not p_t_ratio:

        if M:
            p_t_ratio = (1 + (gamma-1) / 2 * M**2)**(gamma/(gamma-1))

        elif p and p_t:
            p_t_ratio = p_t/p

        elif rho_t_ratio:
            p_t_ratio = rho_t_ratio**(gamma)

        elif T_t_ratio:
            p_t_ratio = T_t_ratio**(gamma/(gamma-1))

    print(f'p_t/p = {p_t_ratio}')
    return p_t_ratio

def get_static_temperature(gamma=1.4,R=287,M=None,u=None,p=None,rho=None,p_t=None,rho_t=None,T_t=None,p_t_ratio=None,rho_t_ratio=None,T_t_ratio=None):

    T = None

    while not T:

        if p and rho:
            T = p/(rho*R)

        elif T_t and M:
            T = T_t/((1 + (gamma-1) / 2 * M**2))

        elif T_t_ratio and T_t:
            T = T_t/T_t_ratio

        elif p_t_ratio and T_t:
            T = T_t/p_t_ratio**((gamma-1)/gamma)

        elif p and p_t and T_t:
            T = T_t/(p_t/p)**((gamma-1)/gamma)

        elif rho_t_ratio and T_t:
            T = T_t/rho_t_ratio**(1/(gamma-1))

        elif rho and rho_t and T_t:
            T = T_t/(rho_t/rho)**(1/(gamma-1))

        elif T_t and (u == 0):
            T = T_t

    print(f'Static temperature = {T}')
    return T

def get_total_temperature(gamma=1.4,R=287,M=None,u=None,p=None,rho=None,T=None,p_t=None,rho_t=None,p_t_ratio=None,rho_t_ratio=None,T_t_ratio=None):

    T_t = None

    while not T_t:

        if p_t and rho_t:
            T_t = p_t/(rho_t*R)

        elif T and M:
            T_t = T*((1 + (gamma-1) / 2 * M**2))

        elif T and T_t:
            T_t = T*T_t_ratio

        elif p_t_ratio and T:
            T_t = T*p_t_ratio**((gamma-1)/gamma)

        elif p and p_t and T:
            T_t = T*(p_t/p)**((gamma-1)/gamma)            

        elif rho_t_ratio and T:
            T_t = T*rho_t_ratio**(1/(gamma-1))

        elif rho and rho_t and T:
            T_t = T*(rho_t/rho_t_ratio)**(1/(gamma-1))

        elif T and (u == 0):
            T_t = T

    print(f'Total temperature = {T_t}')
    return T_t

def get_temperature_ratio(gamma=1.4,R=287,M=None,u=None,p=None,rho=None,T=None,p_t=None,rho_t=None,T_t=None,p_t_ratio=None,rho_t_ratio=None):

    T_t_ratio = None

    while not T_t_ratio:

        if M:
            T_t_ratio = (1 + (gamma-1) / 2 * M**2)

        elif T and T_t:
            T_t_ratio = T_t/T

        elif rho_t_ratio:
            T_t_ratio = rho_t_ratio**(gamma-1)

        elif p_t_ratio:
            T_t_ratio = p_t_ratio**((gamma-1)/gamma)

    print(f'T_t/T = {T_t_ratio}')
    return T_t_ratio

def get_static_density(gamma=1.4,R=287,M=None,u=None,p=None,T=None,p_t=None,rho_t=None,T_t=None,p_t_ratio=None,rho_t_ratio=None,T_t_ratio=None):

    rho = None
    
    while not rho:

        if p and T:
            rho = p/(R*T)

        elif rho_t and M:
            rho = rho_t/((1 + (gamma-1) / 2 * M**2))**(1/(gamma-1))

        elif rho_t and rho_t_ratio:
            rho = rho_t/rho_t_ratio
        
        elif p_t_ratio and rho_t:
            rho = rho_t/p_t_ratio**(1/gamma)

        elif p and p_t and rho_t:
            rho = rho_t/(p_t/p)**(1/gamma)

        elif T_t_ratio and rho_t:
            rho = rho_t/T_t_ratio**(1/(gamma-1))

        elif T and T_t and rho_t:
            rho = rho_t/(T_t/T)**(1/(gamma-1))

        elif rho_t and (u == 0):
            rho = rho_t

    print(f'Static density = {rho}')
    return rho

def get_total_density(gamma=1.4,R=287,M=None,u=None,p=None,rho=None,T=None,p_t=None,T_t=None,p_t_ratio=None,rho_t_ratio=None,T_t_ratio=None):

    rho_t = None
    
    while not rho_t:

        if p_t and T_t:
            rho_t = p_t/(R*T_t)

        elif rho and M:
            rho_t = rho*((1 + (gamma-1) / 2 * M**2))**(1/(gamma-1))

        elif rho and rho_t_ratio:
            rho_t = rho*rho_t_ratio
        
        elif p_t_ratio and rho:
            rho_t = rho*p_t_ratio**(1/gamma)

        elif p and p_t and rho:
            rho_t = rho*(p_t/p)**(1/gamma)
        
        elif T_t_ratio and rho:
            rho_t = rho*T_t_ratio**(1/(gamma-1))

        elif T and T_t and rho:
            rho_t = rho*(T_t/T)**(1/(gamma-1))

        elif rho and (u == 0):
            rho_t = rho

    print(f'Total density = {rho_t}')
    return rho_t

def get_density_ratio(gamma=1.4,R=287,M=None,u=None,p=None,rho=None,T=None,p_t=None,rho_t=None,T_t=None,p_t_ratio=None,T_t_ratio=None):

    rho_t_ratio = None

    while not rho_t_ratio:

        if M:
            rho_t_ratio = (1 + (gamma-1) / 2 * M**2)**(1/(gamma-1))

        elif rho and rho_t:
            rho_t_ratio = rho_t/rho

        elif T_t_ratio:
            rho_t_ratio = T_t_ratio**(1/(gamma-1))

        elif p_t_ratio:
            rho_t_ratio = p_t_ratio**(1/gamma)

    print(f'rho_t/rho = {rho_t_ratio}')
    return rho_t_ratio

if __name__=='__main__':
    #u = get_fluid_velocity(M=2,a=100)
    #u1 = get_fluid_velocity(M=2,T=300)
    #u2 = get_fluid_velocity(M=2,p=101300,rho=1.22)
    #m = get_mach_number(a=100,u=300)
    #m1 = get_mach_number(u=300,T=400)
    #m2 = get_mach_number(u=300,p=101300,rho=1.22)
    foo = CompressibleFlow(p_t_ratio=1.5)
    foo.isentropic_state()
    foo.report_outputs()
    bar = CompressibleFlow(M=1.2,p=120000,T=300)
    bar.isentropic_state()
    bar.report_outputs()
    #one = get_mach_number(T_t_ratio=2)
    #two = get_mach_number(p_t_ratio=2)
    #three = get_mach_number(rho_t_ratio=2)
    print('I am done')