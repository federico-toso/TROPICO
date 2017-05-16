function vehicle = load_vehicle(vtype)

switch vtype
    
    case('hyperion')
        vehicle.Sgross = 200;
        vehicle.gtow = 300e3;
        vehicle.m0 = 30e3;               % Mass (kg)
        vehicle.Tangle = 0*pi/180;       % thrust angle offset with respect to x-y body plane (rad)

    case('IAC16-C1')
        vehicle.gtow = 90035;
        vehicle.m0 = 23189;%23189;               % Mass (kg) 
        vehicle.Sgross = 355;            % Gross wing area(m^2)
        vehicle.Tangle = 6*pi/180;       % thrust angle offset with respect to x-y body plane (rad)         
        vehicle.Rnose = 1;             % Nose radius (m)
        vehicle.Epsilon = 0.8;           % Surface emissivity
        
    case('ISTS')
        vehicle.gtow = 100e3;
        vehicle.m0 = 0.3*vehicle.gtow;%23189;               % Mass (kg) 
        vehicle.Sgross = 1.52535e-3*vehicle.gtow;           % Gross wing area(m^2)  %estimated on 100 ton vehicle
        vehicle.Tangle = 0*pi/180;      % thrust angle offset with respect to x-y body plane (rad)         
        vehicle.Rnose = 1;              % Nose radius (m)
        vehicle.Epsilon = 0.8;          % Surface emissivity
        
    case('AIAAHYP17')
        vehicle.gtow = 90035;
        vehicle.m0 = 0.25443*vehicle.gtow;%23189;               % Mass (kg) 
        vehicle.Sgross = 4.05e-3*vehicle.gtow;           % Gross wing area(m^2)  %estimated on 100 ton vehicle
        vehicle.Tangle = 0*pi/180;      % thrust angle offset with respect to x-y body plane (rad)         
        vehicle.Rnose = 1;              % Nose radius (m)
        vehicle.Epsilon = 0.8;          % Surface emissivity
        
    case('IAC16_TUG')
        vehicle.gtow = 1e4;
        vehicle.m0 = 1e2;               % Mass (kg) 
        vehicle.Sgross = 1;            % Gross wing area(m^2)
        vehicle.Tangle = 0*pi/180;       % thrust angle offset with respect to x-y body plane (rad)         
        vehicle.Rnose = 0.5;             % Nose radius (m)
        vehicle.Epsilon = 0.6;           % Surface emissivity       
        
    case('UPPER_STAGE')
        vehicle.gtow = 1e3;
        vehicle.m0 = 1e2;               % Mass (kg) 
        vehicle.Sgross = 1;            % Gross wing area(m^2)
        vehicle.Tangle = 0*pi/180;       % thrust angle offset with respect to x-y body plane (rad)         
        vehicle.Rnose = 0.5;             % Nose radius (m)
        vehicle.Epsilon = 0.6;           % Surface emissivity
        
    case('test')
        vehicle.gtow = 100e3;
        vehicle.m0 = 1e3;               % Mass (kg) 
        vehicle.Sgross = 355;            % Gross wing area(m^2)
        vehicle.Tangle = 0*pi/180;       % thrust angle offset with respect to x-y body plane (rad)         
        vehicle.Rnose = 0.5;             % Nose radius (m)
        vehicle.Epsilon = 0.6;           % Surface emissivity
        
    case('falcon9_S1+2')
        vehicle.gtow = 550e3;
        vehicle.m0 = 34e3;               % Mass (kg) 
        vehicle.Sgross = 3*3*pi;            % Gross wing area(m^2)
        vehicle.Tangle = 0;       % thrust angle offset with respect to x-y body plane (rad)         
        vehicle.Rnose = 2;             % Nose radius (m)
        vehicle.Epsilon = 0.6;           % Surface emissivity
        
    case('falcon9_S2')
        vehicle.gtow = 125e3;
        vehicle.m0 = 4e3;               % Mass (kg) 
        vehicle.Sgross = 2*2*pi;            % Gross wing area(m^2)
        vehicle.Tangle = 0;       % thrust angle offset with respect to x-y body plane (rad)         
        vehicle.Rnose = 2;             % Nose radius (m)
        vehicle.Epsilon = 0.6;           % Surface emissivity
  
end