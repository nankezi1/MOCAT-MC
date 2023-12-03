function [objclass_out] = objclass2int(objclass_in,direction)
%OBJCLASS2INT converts objclass_in string into an integer (direction=1),
%or vice versa (direction=2)
%         'Payload' = 1
%         'Payload Mission Related Object' = 2
%         'Payload Fragmentation Debris' = 3
%         'Payload Debris' = 4
%         'Rocket Body' = 5
%         'Rocket Mission Related Object' = 6
%         'Rocket Fragmentation Debris' = 7
%         'Rocket Debris' = 8
%         'Debris' (new) = 9
%         'Other Debris' = 10
%         'Unknown' or [] = 11
%         'Untracked Debris' = 12

possible_objclass = {'Payload';...
                     'Payload Mission Related Object';...
                     'Payload Fragmentation Debris';...
                     'Payload Debris';...
                     'Rocket Body';...
                     'Rocket Mission Related Object';...
                     'Rocket Fragmentation Debris';...
                     'Rocket Debris';...
                     'Debris';...
                     'Other Debris';...
                     'Unknown';...
                     'Untracked Debris'};

if direction==1 %string to integer
    
    if isequal(objclass_in,possible_objclass{1})
        objclass_out = 1;
    elseif isequal(objclass_in,possible_objclass{2})
        objclass_out = 2;
    elseif isequal(objclass_in,possible_objclass{3})
        objclass_out = 3;
    elseif isequal(objclass_in,possible_objclass{4})
        objclass_out = 4;
    elseif isequal(objclass_in,possible_objclass{5})
        objclass_out = 5;
    elseif isequal(objclass_in,possible_objclass{6})
        objclass_out = 6;
    elseif isequal(objclass_in,possible_objclass{7})
        objclass_out = 7;
    elseif isequal(objclass_in,possible_objclass{8})
        objclass_out = 8;
    elseif isequal(objclass_in,possible_objclass{9})
        objclass_out = 9;
    elseif isequal(objclass_in,possible_objclass{10})
        objclass_out = 10;
    elseif isequal(objclass_in,possible_objclass{11}) || isempty(objclass_in)
        objclass_out = 11;
    elseif isequal(objclass_in,possible_objclass{12})
        objclass_out = 12;
    else 
        error('objclass_in ''%s'' did not match any of the pre-determined options. Please review this function.',objclass_in)
    end
    
elseif direction==2 %integer to string
    if length(objclass_in)==1
        objclass_out = possible_objclass{objclass_in}; %output is a string
    else
        objclass_out = possible_objclass(objclass_in); %output is a cell array [length(objclass_in) x 1]
    end
else
    error('direction does not match 1 or 2')
end

