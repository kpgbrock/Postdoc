function [lenStr] = getLengthStringWithoutCharacter(inStr)

charToRemove = '_';
inStr(inStr==charToRemove) = '';
lenStr = length(inStr);

end