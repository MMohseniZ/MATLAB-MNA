clc; clear; close all;

syms t s
assume(t, 'real'); % this assumes the t as real value as it is
prompt = 'Enter netList file name : ';
inputFileName = input(prompt, "s"); % get input file name from command window
clc;
prompt = 'Enter output file name : ';
outputFileName = input(prompt, "s"); % get output file name from command window
clc;
fileID = fopen(inputFileName);
compLine = fgetl(fileID);
digits(2);
elements = {};

Nodes = [];

while ischar(compLine) % this while loop read the elements and copies them to an array with element type which defined as class
    currentComponent = split(compLine, ',');
    
    if ~ (any(Nodes(:) == str2double(currentComponent(3))))
        Nodes(end + 1) = str2double(currentComponent(3));
    end

    if ~ (any(Nodes(:) == str2double(currentComponent(4))))
        Nodes(end + 1) = str2double(currentComponent(4));
    end
    if strcmp(currentComponent(1), 'ML')
        if ~(any(Nodes(:) == str2double(currentComponent(6))))
            Nodes(end + 1) = str2double(currentComponent(6));
        end
    
        if ~(any(Nodes(:) == str2double(currentComponent(7))))
            Nodes(end + 1) = str2double(currentComponent(7));
        end
    end
    
    if (strcmp(currentComponent(1), 'R') || strcmp(currentComponent(1), 'C') || strcmp(currentComponent(1), 'L') || strcmp(currentComponent(1), 'V') || strcmp(currentComponent(1), 'I'))
        if strcmp(currentComponent(1), 'V') || strcmp(currentComponent(1), 'I')
            val = laplace(str2sym(currentComponent(5)), s);
            IC = [];
        elseif strcmp(currentComponent(1), 'R')
            val = 1./str2double(currentComponent(5));
            IC = 0;
        elseif strcmp(currentComponent(1), 'C')
            val = str2double(currentComponent(5)) .* s;
            IC = str2double(currentComponent(5)) * str2double(currentComponent(6));
        elseif strcmp(currentComponent(1), 'L')
            val = str2double(currentComponent(5)) .* s;
            IC = str2double(currentComponent(5)) * str2double(currentComponent(6));
        end
        elements{end + 1} = element(currentComponent(1), currentComponent(2), [str2double(currentComponent(3)), str2double(currentComponent(4))], [], val, IC);
    end
    if strcmp(currentComponent(1), 'ML')
        elements{end + 1} = element('ML', currentComponent(2), [str2double(currentComponent(3)), str2double(currentComponent(4)), str2double(currentComponent(6)), str2double(currentComponent(7))], [], [str2double(currentComponent(5)), str2double(currentComponent(8)), str2double(currentComponent(9))], []);
    end
    if (strcmp(currentComponent(1), 'Z') || strcmp(currentComponent(1), 'H') || strcmp(currentComponent(1), 'Y') || strcmp(currentComponent(1), 'T'))
        gain = str2sym(currentComponent(7));
        elements{end + 1} = element(currentComponent(1) , currentComponent(2), [str2double(currentComponent(3)), str2double(currentComponent(4))], [str2double(currentComponent(5)), str2double(currentComponent(6))], gain, []);
    end
    compLine = fgetl(fileID);
end
fclose(fileID);
numOfNodes = max(Nodes, [], 'all');

A = sym(zeros(numOfNodes)); % create the left hand side matrix of equations
B = sym(zeros(numOfNodes, 1)); % create the right hand side matrix of equations

for n = 1:length(elements) % complete the A matrix regarding to elements which there is no need to know their current directly
    type = elements{n}.type;
    if strcmp(type, 'R') || strcmp(type, 'C')
        i = elements{n}.nodes(1);
        j = elements{n}.nodes(2);
        val = elements{n}.value;
        IC = elements{n}.IC;
        if i == 0
            A(j, j) = A(j, j) + val;
            B(j) = B(j) - IC;
        elseif j == 0
            A(i, i) = A(i, i) + val;
            B(i) = B(i) + IC;
        else
            A(i, i) = A(i, i) + val;
            A(j, j) = A(j, j) + val;
            A(i, j) = A(i, j) - val;
            A(j, i) = A(j, i) - val;
            B(j) = B(j) - IC;
            B(i) = B(i) + IC;
        end
    elseif strcmp(type, 'I')
        i = elements{n}.nodes(1);
        j = elements{n}.nodes(2);
        val = elements{n}.value;
        if i == 0
            B(j) = B(j) + val;
        elseif j == 0
            A(i, i) = A(i, i) + val;
            B(i) = B(i) - val;
        else
            B(j) = B(j) + val;
            B(i) = B(i) - val;
        end
    elseif strcmp(type, 'H')
        k = elements{n}.nodes(1);
        l = elements{n}.nodes(2);
        i = elements{n}.snodes(1);
        j = elements{n}.snodes(2);
        val = elements{n}.value;
        if i == 0
            if k == 0
                A(l, j) = A(l, j) + val;
            elseif l == 0
                A(k, j) = A(k, j) - val;
            else
                A(l, j) = A(l, j) + val;
                A(k, j) = A(k, j) - val;
            end
        elseif j == 0
            if k == 0
                A(l, i) = A(l, i) - val;
            elseif l == 0
                A(k, i) = A(k, i) + val;
            else
                A(l, i) = A(l, i) - val;
                A(k, i) = A(k, i) + val;
            end
        else
            if k == 0
                A(l, i) = A(l, i) - val;
                A(l, j) = A(l, j) + val;
            elseif l == 0
                A(k, i) = A(k, i) + val;
                A(k, j) = A(k, j) - val;
            else
                A(l, i) = A(l, i) - val;
                A(k, i) = A(k, i) + val;
                A(l, j) = A(l, j) + val;
                A(k, j) = A(k, j) - val;
            end
        end
    end
end

for n = 1:length(elements) % complete the A matrix regarding to elements which we have to know their current directly to solve the node equations
    type = elements{n}.type;
    if strcmp(type, 'L')
        i = elements{n}.nodes(1);
        j = elements{n}.nodes(2);
        val = elements{n}.value;
        IC = elements{n}.IC;
        A = [A, sym(zeros(length(A(:,1)), 1));
            sym(zeros(1, length(A(:,1)))), sym(0)];
        B = [B; sym(0)];
        A(end, end) = -val;
        B(end) = -IC;
        if i == 0
            A(end, j) = -1 + A(end, j);
            A(j, end) = -1 + A(j, end);
        elseif j == 0
            A(end, i) = 1 + A(end, i);
            A(i, end) = 1 + A(i, end);
        else
            A(end, j) = -1 + A(end, j);
            A(j, end) = -1 + A(j, end);
            A(end, i) = 1 + A(end, i);
            A(i, end) = 1 + A(i, end);
        end
    elseif strcmp(type, 'V')
        i = elements{n}.nodes(1);
        j = elements{n}.nodes(2);
        val = elements{n}.value;
        IC = elements{n}.IC;
        A = [A, sym(zeros(length(A(:,1)), 1));
            sym(zeros(1, length(A(:,1)))), sym(0)];
        B = [B; val];
        if i == 0
            A(end, j) = -1 + A(end, j);
            A(j, end) = -1 + A(j, end);
        elseif j == 0
            A(end, i) = 1 + A(end, i);
            A(i, end) = 1 + A(i, end);
        else
            A(end, j) = -1 + A(end, j);
            A(j, end) = -1 + A(j, end);
            A(end, i) = 1 + A(end, i);
            A(i, end) = 1 + A(i, end);
        end
    elseif strcmp(type, 'ML')
        i = elements{n}.nodes(1);
        j = elements{n}.nodes(2);
        k = elements{n}.nodes(3);
        l = elements{n}.nodes(4);
        val = elements{n}.value;
        A = [A, sym(zeros(length(A(:,1)), 2));
            sym(zeros(2, length(A(:,1)))), sym(zeros(2))];
        B = [B; sym(zeros(2, 1))];
        A(end - 1, end - 1) = -elements{n}.value(1)*s;
        A(end - 1, end) = -elements{n}.value(3)*s;
        A(end, end - 1) = -elements{n}.value(3)*s;
        A(end, end) = -elements{n}.value(2)*s;
        if i == 0
            A(j, end - 1) = -1 + A(j, end - 1);
            A(end - 1, j) = -1 + A(end - 1, j);
        elseif j == 0
            A(i, end - 1) = 1 + A(i, end - 1);
            A(end - 1, i) = 1 + A(end - 1, i);
        else
            A(j, end - 1) = -1 + A(j, end - 1);
            A(end - 1, j) = -1 + A(end - 1, j);
            A(i, end - 1) = 1 + A(i, end - 1);
            A(end - 1, i) = 1 + A(end - 1, i);
        end

        if k == 0
            A(l, end) = -1 + A(l, end);
            A(end, l) = -1 + A(end, l);
        elseif l == 0
            A(k, end) = 1 + A(k, end);
            A(end, k) = 1 + A(end, k);
        else
            A(l, end) = -1 + A(l, end);
            A(end, l) = -1 + A(end, l);
            A(k, end) = 1 + A(k, end);
            A(end, k) = 1 + A(end, k);
        end
    elseif strcmp(type, 'Z') || strcmp(type, 'T')
        k = elements{n}.nodes(1);
        l = elements{n}.nodes(2);
        i = elements{n}.snodes(1);
        j = elements{n}.snodes(2);
        val = elements{n}.value;
        A = [A, sym(zeros(length(A(:,1)), 1));
            sym(zeros(1, length(A(:,1)))), sym(0)];
        B = [B; sym(0)];
        if strcmp(type, 'Z')
            if i == 0
                A(end, j) = val + A(end, j);
            elseif j == 0
                A(end, i) = -val + A(end, i);
            else
                A(end, j) = val + A(end, j);
                A(end, i) = -val + A(end, i);
            end
    
            if k == 0
                A(l, end) = -1 + A(l, end);
                A(end, l) = -1 + A(end, l);
            elseif l == 0
                A(end, k) = 1 + A(end, k);
                A(k, end) = 1 + A(k, end);
            else
                A(l, end) = -1 + A(l, end);
                A(end, l) = -1 + A(end, l);
                A(end, k) = 1 + A(end, k);
                A(k, end) = 1 + A(k, end);
            end
        else
            if i == 0
                A(j, end) = -1 + A(j, end);
                A(end, j) = -1 + A(end, j);
            elseif j == 0
                A(end, i) = 1 + A(end, i);
                A(i, end) = 1 + A(i, end);
            else
                A(j, end) = -1 + A(j, end);
                A(end, j) = -1 + A(end, j);
                A(end, i) = 1 + A(end, i);
                A(i, end) = 1 + A(i, end);
            end
    
            if k == 0
                A(l, end) = -val + A(l, end);
            elseif l == 0
                A(k, end) = val + A(k, end);
            else
                A(l, end) = -val + A(l, end);
                A(k, end) = val + A(k, end);
            end
        end

    elseif strcmp(type, 'Y')
        k = elements{n}.nodes(1);
        l = elements{n}.nodes(2);
        i = elements{n}.snodes(1);
        j = elements{n}.snodes(2);
        val = elements{n}.value;
        A = [A, sym(zeros(length(A(:,1)), 2));
            sym(zeros(2, length(A(:,1)))), sym(zeros(2))];
        B = [B; sym(zeros(2, 1))];
        A(end, end - 1) = -val;
        if i == 0
            A(j, end - 1) = -1 + A(j, end - 1);
            A(end - 1, j) = -1 + A(end - 1, j);
        elseif j == 0
            A(i, end - 1) = 1 + A(i, end - 1);
            A(end - 1, i) = 1 + A(end - 1, i);
        else
            A(j, end - 1) = -1 + A(j, end - 1);
            A(end - 1, j) = -1 + A(end - 1, j);
            A(i, end - 1) = 1 + A(i, end - 1);
            A(end - 1, i) = 1 + A(end - 1, i);
        end

        if k == 0
            A(l, end) = -1 + A(l, end);
            A(end, l) = -1 + A(end, l);
        elseif l == 0
            A(k, end) = 1 + A(k, end);
            A(end, k) = 1 + A(end, k);
        else
            A(l, end) = -1 + A(l, end);
            A(end, l) = -1 + A(end, l);
            A(k, end) = 1 + A(k, end);
            A(end, k) = 1 + A(end, k);
        end
    end
end

X = linsolve(A, B); % solve the system of equations
X = simplify(simplifyFraction(X)); % simplify the expressions of results
X = [sym(0); X]; % enter the ground node to the matrix
IDN = 1; % number of current dependent elements
numOfNodes = numOfNodes + 1;

for n = 1:length(elements) % in this for we assign the voltage, current and power of each element in elements array
    type = elements{n}.type;
    if ~ strcmp(type, 'ML')
        elements{n}.voltage = X(elements{n}.nodes(1) + 1) - X(elements{n}.nodes(2) + 1);
        if strcmp(type, 'C') || strcmp(type, 'R')
            elements{n}.current = elements{n}.value * elements{n}.voltage;
        elseif strcmp(type, 'L') || strcmp(type, 'V') || strcmp(type, 'Z') || strcmp(type, 'T')
            elements{n}.current = X(numOfNodes + IDN);
            IDN = IDN + 1;
        elseif strcmp(type, 'I')
            elements{n}.current = elements{n}.value;
        elseif strcmp(type, 'H')
            elements{n}.current = elements{n}.value * (X(elements{n}.snodes(1) + 1) - X(elements{n}.snodes(2) + 1));
        elseif strcmp(type, 'Y')
            elements{n}.current = X(numOfNodes + IDN + 1);
            IDN = IDN + 2;
        end
        elements{n}.voltage = simplify(real(ilaplace(elements{n}.voltage)));
        elements{n}.current = simplify(real(ilaplace(elements{n}.current)));
        elements{n}.power = simplify(elements{n}.voltage * elements{n}.current);
    else
        elements{n}.voltage = sym(zeros(1,2));
        elements{n}.current = sym(zeros(1,2));
        elements{n}.voltage(1) = simplify(real(ilaplace(X(elements{n}.nodes(1) + 1) - X(elements{n}.nodes(2) + 1))));
        elements{n}.voltage(2) = simplify(real(ilaplace(X(elements{n}.nodes(3) + 1) - X(elements{n}.nodes(4) + 1))));
        elements{n}.current(1) = simplify(real(ilaplace(X(numOfNodes + IDN))));
        elements{n}.current(2) = simplify(real(ilaplace(X(numOfNodes + IDN + 1))));
        elements{n}.power = simplify(elements{n}.voltage(1)*elements{n}.current(1) + elements{n}.voltage(2)*elements{n}.current(2));
        IDN = IDN + 2;
    end
end

fileID = fopen(outputFileName,'w');

for n = 1:length(elements) % in this for we generate the output text file
    type =  elements{n}.type;
    if strcmp(type, 'ML')
        name = elements{n}.name;
        v1 = real(vpa(elements{n}.voltage(1)));
        v2 = real(vpa(elements{n}.voltage(2)));
        i1 = real(vpa(elements{n}.current(1)));
        i2 = real(vpa(elements{n}.current(2)));
        p = real(vpa(elements{n}.power));
        fprintf(fileID,'<%s><%s><%s><%s><%s><%s>\n', [name v1 i1 v2 i2 p]);
    else
        name = elements{n}.name;
        v = real(vpa(elements{n}.voltage));
        i = real(vpa(elements{n}.current));
        p = real(vpa(elements{n}.power));
        fprintf(fileID,'<%s><%s><%s><%s>\n', [name v i p]);
    end
end
fclose(fileID);

disp('output file generated ...');
in = ' ';
while ~strcmp(in, 'end') % this is the while loop which is used to get input repeatedly from user to plot informations of elements
    in = input('enter name of the element or stop simulation using "end" : ', "s");
    clc;
    b = 0;
    for i = 1:length(elements)
        if strcmp(elements{i}.name, in)
            interval = str2num(input('Enter the time interval (ex : "0 3") : ', "s"));
            dt = (interval(2) - interval(1))/1000;
            interval = [interval(1) + dt, interval(2) - dt];
            e = elements{i};
            b = 1;
            if strcmp(e.type, 'ML')
                figure
                subplot(3, 2, 1)
                fplot(e.voltage(1), interval, 'b', 'LineWidth', 2);
                xlabel('Time(s)', 'Interpreter', 'latex');
                ylabel('V1(Volts)', 'Interpreter', 'latex');
                title(strcat('V1', '(', e.name, ')'), 'Interpreter', 'latex');
                grid minor;
                subplot(3, 2, 2)
                fplot(e.current(1), interval, 'r', 'LineWidth', 2);
                xlabel('Time(s)', 'Interpreter', 'latex');
                ylabel('I1(Amps)', 'Interpreter', 'latex');
                title(strcat('I1', '(', e.name, ')'), 'Interpreter', 'latex');
                grid minor;
                subplot(3, 2, 3)
                fplot(e.voltage(2), interval, 'b', 'LineWidth', 2);
                xlabel('Time(s)', 'Interpreter', 'latex');
                ylabel('V2(Volts)', 'Interpreter', 'latex');
                title(strcat('V2', '(', e.name, ')'), 'Interpreter', 'latex');
                grid minor;
                subplot(3, 2, 4)
                fplot(e.current(2), interval, 'r', 'LineWidth', 2);
                xlabel('Time(s)', 'Interpreter', 'latex');
                ylabel('I2(Amps)', 'Interpreter', 'latex');
                title(strcat('I2', '(', e.name, ')'), 'Interpreter', 'latex');
                grid minor;
                subplot(3, 2, [5 6])
                fplot(e.power, interval, 'green', 'LineWidth', 2);
                xlabel('Time(s)', 'Interpreter', 'latex');
                ylabel('P(Watts)', 'Interpreter', 'latex');
                title(strcat('P', '(', e.name, ')'), 'Interpreter', 'latex');
                grid minor;
            else
                subplot(3, 1, 1)
                fplot(e.voltage, interval, 'b', 'LineWidth', 2);
                xlabel('Time(s)', 'Interpreter', 'latex');
                ylabel('V(Volts)', 'Interpreter', 'latex');
                title(strcat('V', '(', e.name, ')') , 'Interpreter', 'latex');
                grid minor;
                subplot(3, 1, 2)
                fplot(e.current, interval, 'r', 'LineWidth', 2);
                xlabel('Time(s)', 'Interpreter', 'latex');
                ylabel('I(Amps)', 'Interpreter', 'latex');
                title(strcat('I', '(', e.name, ')'), 'Interpreter', 'latex');
                grid minor;
                subplot(3, 1, 3)
                fplot(e.power, interval, 'green', 'LineWidth', 2);
                xlabel('Time(s)', 'Interpreter', 'latex');
                ylabel('P(Watts)', 'Interpreter', 'latex');
                title(strcat('P', '(', e.name, ')'), 'Interpreter', 'latex');
                grid minor;
            end
        end
    end
    if b == 0 && ~strcmp(in, 'end')
        disp('there is no element with this name in circuit');
    end
end

clc; close all;
