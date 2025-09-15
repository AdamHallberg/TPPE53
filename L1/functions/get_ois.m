function [OISs, option_data] = real_excel_file(filename)
    ois_data = readtable(filename);
    d = ois_data(2:end, [5 6]);
    % Read maturities (col D) and mid rates (col F) from Excel
    maturities = floor(d(:,"Var5")-d(1,"Var5"))./365; % [years]
    mids       = d(:,2).*0.01;
    n = max(size(maturities));

    OISs = [maturities, mids];


    option_data = readtable("ois_data.xlsx", "Sheet", "filtered_option_data");
    option_data = table2array(option_data);

    % 4 is "arbitrary"
    option_data(:, 3) = option_data(:, 3) - option_data(1, 3) + 4; % get time diff
    option_data(:, 3) = option_data(:, 3) / 365; % conver to fraction of years
    
    
end
