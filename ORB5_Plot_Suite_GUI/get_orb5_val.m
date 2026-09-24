function value = get_orb5_val(blockName,paramName,filename)
%GET_ORB5_VAL Parse simple namelists without evaluating numeric expressions.
if nargin < 3, filename = fullfile(pwd,'input'); end
content = fileread(filename);
[content,strings] = protectStrings(content);
block = regexp(content,['(?is)&' regexptranslate('escape',blockName) ...
    '(?!\w)(.*?)(?:/|&end)'],'tokens','once');
if isempty(block), error('ORB5:MissingBlock','Block &%s not found in %s.',blockName,filename); end
token = regexp(block{1},['(?is)(?<!\w)' regexptranslate('escape',paramName) ...
    '\s*=\s*(.*?)(?=,?\s*(?<!\w)[A-Za-z]\w*\s*=|$)'],'tokens','once');
if isempty(token), error('ORB5:MissingParameter','Parameter %s missing from &%s.',paramName,blockName); end
raw = regexprep(strtrim(token{1}),',\s*$','');
if strcmpi(raw,'.true.'), value = true;
elseif strcmpi(raw,'.false.'), value = false;
elseif startsWith(raw,'ORB5STRINGTOKEN')
    index = sscanf(raw,'ORB5STRINGTOKEN%d');
    value = strings{index};
else
    parts = regexp(strtrim(strrep(raw,',',' ')),'\s+','split');
    value = str2double(regexprep(parts,'[dD]','e'));
    if any(isnan(value)), error('ORB5:UnsupportedValue','Unsupported value for %s: %s',paramName,raw); end
end
end

function [out,strings]=protectStrings(content)
% Quoted ! and / are data, not comment/block delimiters.
out=''; strings={}; k=1;
while k<=numel(content)
    ch=content(k);
    if ch=='!'
        while k<=numel(content) && ~ismember(content(k),[char(10),char(13)]), k=k+1; end
    elseif any(ch==['''' '"'])
        quote=ch; k=k+1; value=''; closed=false;
        while k<=numel(content)
            if content(k)==quote
                if k<numel(content) && content(k+1)==quote
                    value(end+1)=quote; k=k+2; continue;
                end
                k=k+1; closed=true; break;
            end
            value(end+1)=content(k); k=k+1;
        end
        if ~closed, error('ORB5:Namelist','Unterminated quoted string.'); end
        strings{end+1}=value;
        out=[out sprintf('ORB5STRINGTOKEN%d',numel(strings))]; %#ok<AGROW>
    else
        out(end+1)=ch; k=k+1;
    end
end
end

