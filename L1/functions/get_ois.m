function OISs = get_ois(filename)

    % Read maturities (col D) and mid rates (col F) from Excel
    maturities = readmatrix(filename, 'Range', 'D3:D24');
    mids       = readmatrix(filename, 'Range', 'F3:F24');

    % Preallocate struct array
    n = length(maturities);
    OISs = struct('Maturity', cell(1, n), 'Mid', cell(1, n));

    % Fill struct array
    for i = 1:n
        OISs(i).Maturity = maturities(i);
        OISs(i).Mid      = mids(i);
    end
end