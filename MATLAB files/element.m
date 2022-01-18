classdef element
    properties
        type % type of element
        name % name of element
        nodes % nodes of element
        snodes % source nodes for VCVS, CCVS, VCCS, CCCS
        value % values needed for element charactrization
        IC % initial condition
        voltage % voltage of element
        current % current of element
        power % the power of element
    end

    methods
        function obj = element(type, name, nodes, snodes, value, IC)
            obj.type = type;
            obj.name = name;
            obj.nodes = nodes;
            obj.snodes = snodes;
            obj.value = value;
            obj.IC = IC;
        end
    end
end
