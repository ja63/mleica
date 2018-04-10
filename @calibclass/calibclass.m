classdef calibclass < handle
    %     clean all
    %     clear all
    %     load labdata.mat
    %     mo = myclass;
    %     load('labdata.mat')
    %     mo.calibrate(A,B)
    %     mo.lc2rc(B)
    
    properties (Constant)
        fs = 100;
        jaDebug = true;
        % Hidden functions
        jl = @(x,c) line(x(:,1),x(:,2),x(:,3),'color',c);
        jdist = @(q) sqrt(sum(diff(q).^2,2));
        js = @(x,c) scatter3(x(:,1),x(:,2),x(:,3),c);
        jbp = @(x) calibclass.ballplot(x(:,1),x(:,2),x(:,3),'r',20);
        jbpcr = @(x,c,r) calibclass.ballplot(x(:,1),x(:,2),x(:,3),c,r);
        cl = @() camlight; % Usage cl()
        fw = @() view([1 0 0]);
        lw = @() view([0 -1 0]);
        rw = @() view([0 1 0]);
        tw = @() view([0 0 1]);
        jv = @(q,fs) medfilt1(sqrt(sum(diff(q).^2,2))*fs,5);
    end
    properties
        % Variables for calculus
        RM
        TMA
        TMB
    end
    methods
        function obj1 = myclass(varargin)
            if nargin == 0
                
            end
        end
        function calibrate(this,varargin)
            % A is robot coordinates
            % B is Leica coordinates
            if nargin == 1
                pth = getpref('calibclass','prgpath',[cd '\']);
                [fn,pth]=uigetfile({'*.mod;*.prg','Rapid files (*.mod, *.prg)'},'Get rapid file',pth);
                switch class(fn)
                    case 'char'
                        disp('Hej');
                        setpref('calibclass','prgpath',pth)
                    case 'double'
                        disp('User cancelled');
                        return
                end
                A = this.parserapid([pth,fn],10:10:100);
                %                 [fn,pth]=uigetfile('*.out','Get leica file');
                pth = getpref('calibclass','outpath',[cd '\']);
                [fn,pth]=uigetfile([pth '*.out'],'Get raw data file');
                switch class(fn)
                    case 'char'
                        disp('Hej');
                        setpref('calibclass','outpath',pth)
                    case 'double'
                        disp('User cancelled');
                        return
                end
                raw = this.jload([pth,fn]);
                [p,~,~] = this.fndzvp(raw);
                B = permute(mean(p),[3,2,1]);
            elseif nargin == 3
                A = this.parserapid(varargin{1},10:10:100);
                raw = this.jload(varargin{2});
                [p,~,~] = this.fndzvp(raw);
                B = permute(mean(p),[3,2,1]);
            else
                return
            end
            this.TMA = repmat(mean(A),size(A,1),1); % traslation to origin
            origA = A - this.TMA;
            this.TMB = repmat(mean(B),size(B,1),1);
            origB = B-this.TMB;
            R = origB\origA;
            [u,~,v] = svd(R);
            this.RM = u*eye(3)*v';
            
            % Jens debuging
            if this.jaDebug
                ortoA = origB*this.RM;
                AAA = ortoA + repmat(mean(A),size(A,1),1);
                figure(6)
                clf(6)
                
                this.jl(A,'r');xlabel('X (mm)');ylabel('Y (mm)');zlabel('Z (mm)');
                hold on
                this.jl(AAA,'b');
                grid;axis equal;
                hold off
            end
            % End jens debugging
        end
        function Rc = lc2rc(this,Lc) %transformation from leica coordinates to robot coordinates
            
            TLc = Lc-repmat(this.TMB(1,:),size(Lc,1),1);
            RLc = TLc*this.RM;
            Rc = RLc+repmat(this.TMA(1,:),size(Lc,1),1);
            
            % Jens debugging
            if false
                figure;hold on;
                this.jl(Rc,'c');
                xlabel('X');
                ylabel('Y');
                zlabel('Z');
            end
            % End Jens debugging
        end
    end
    methods (Static , Access = private)
        function jvr(in,cp1,cp2,tail)
            
            theLine = [cp1;cp2];
            
            % Calculate the normal vector for line P2 -> P1
            nP2P1 = diff(theLine);
            enP2P1 = nP2P1/norm(nP2P1);
            
            % Rotate from a to b
            a = enP2P1;     % unit vectors
            b = [0 0 1];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % see : http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
            ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
            RU = @(A,B) eye(3) + ssc(cross(A,B)) + ...
                ssc(cross(A,B))^2*(1-dot(A,B))/(norm(cross(A,B))^2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            R = RU(a,b);
            
            
            %jl = @(x) line(x(:,1),x(:,2),x(:,3),'color','b');
            jl = @(x,c) line(x(:,1),x(:,2),x(:,3),'color',c);
            % jl(in,'b')
            % hold on
            % jl([R*in']','g')
            
            tjo = @(v) [R*v']';
            
            theLineNew = tjo(theLine);  % Ideal straight line ,
            
            % Rotated to align with z-axis
            newXYZ = cellfun(tjo,in,'UniformOutput',false);
            % Offset mean x,y,z to origo
            %jopp = cellfun(@mean,newXYZ,'UniformOutput',false)
            m = mean(theLineNew);       % Gets centre point
            figure;hold on;
            for i = 1:length(newXYZ)
                newXYZ{i} = newXYZ{i} - repmat(m,size(newXYZ{i},1),1);
                % Mod
                qq = newXYZ{i};
                z = qq(:,3);
                mi = min(z);
                ma = max(z);
                olol = (z < mi+tail) | (z > ma-tail);
                ids = ones(1,size(qq,1));
                ids([find(~olol,1,'first'):find(~olol,1,'last')]) = 0;
                plot3(qq(ids==1,1),qq(ids==1,2),qq(ids==1,3),'r.');
                plot3(qq(ids==0,1),qq(ids==0,2),qq(ids==0,3),'b');
                % End Mod
                % jl(newXYZ{i},'b');
            end
            jl(theLineNew-repmat(m,2,1),'g');
            hold off
        end
        function [I,check]=plane_line_intersect(n,V0,P0,P1)
            %plane_line_intersect computes the intersection of a plane and a segment(or
            %a straight line)
            % Inputs:
            %       n: normal vector of the Plane
            %       V0: any point that belongs to the Plane
            %       P0: end point 1 of the segment P0P1
            %       P1:  end point 2 of the segment P0P1
            %
            %Outputs:
            %      I    is the point of interection
            %     Check is an indicator:
            %      0 => disjoint (no intersection)
            %      1 => the plane intersects P0P1 in the unique point I
            %      2 => the segment lies in the plane
            %      3=>the intersection lies outside the segment P0P1
            %
            % Example:
            % Determine the intersection of following the plane x+y+z+3=0 with the segment P0P1:
            % The plane is represented by the normal vector n=[1 1 1]
            % and an arbitrary point that lies on the plane, ex: V0=[1 1 -5]
            % The segment is represented by the following two points
            % P0=[-5 1 -1]
            %P1=[1 2 3]
            % [I,check]=plane_line_intersect([1 1 1],[1 1 -5],[-5 1 -1],[1 2 3]);
            
            %This function is written by :
            %                             Nassim Khaled
            %                             Wayne State University
            %                             Research Assistant and Phd candidate
            %If you have any comments or face any problems, please feel free to leave
            %your comments and i will try to reply to you as fast as possible.
            
            I=[0 0 0];
            u = P1-P0;
            w = P0 - V0;
            D = dot(n,u);
            N = -dot(n,w);
            check=0;
            if abs(D) < 10^-10        % The segment is parallel to plane
                if N == 0           % The segment lies in plane
                    check=2;
                    return
                else
                    check=0;       %no intersection
                    return
                end
            end
            
            %compute the intersection parameter
            sI = N / D;
            I = P0+ sI.*u;
            
            if (sI < 0 || sI > 1)
                check= 3;          %The intersection point  lies outside the segment, so there is no intersection
            else
                check=1;
            end
            
            
            
        end
        function jpth = jgetfile(varargin)
            if nargin == 2
                pthname = varargin{1};
                ext = varargin{2};
                pth = getpref('calibclass',pthname,[cd '\']);
            elseif nargin == 1
                pthname = 'def';
                ext = varargin{1};
                pth = getpref('calibclass',pthname,[cd '\']);
            else
                pthname = 'def';
                ext = '*';
                pth = getpref('calibclass',pthname,[cd '\']);
            end
            [fn,pth]=uigetfile([pth '*.' ext],'Get rapid file');
            switch class(fn)
                case 'char'
                    disp('Hej');
                    setpref('calibclass',pthname,pth);
                    jpth = [pth,fn];
                case 'double'
                    disp('User cancelled');
                    jpth = 0;
                    return
            end
        end
        function out = parserapid(fn,seq)
            % Fuction returns coordinated within rapid program
            % Input arguments fn = <filename>
            % seq : 10:10:100 gives p10, p20,....,p100
            % Usage :
            %   p = parserapid(['C:\Temp\Koord_synk.prg'], 10:10:100)
            %
            fid = fopen(fn);
            qq = fread(fid,'char');
            fclose(fid);
            qq = char(qq')
            
            p = [];
            n = 1;
            for i = seq
                p(n,:) = str2num(cell2mat(regexp(qq, ...
                    ['(?<=p' num2str(i) ':=[[)[-\d]+\,[-\d]+,[-\d]+'], ...
                    'match')))
                n = n + 1;
            end
            out = p;
        end
    end
    methods (Static)
        function varargout = pointstabtime(qq,theLim)
            % Find out the point stabilisation time according to ISO xxxx and a lot
            % else....  Usage:
            %   [PSt,PSo]=pointstabtime(qq,theLim);
            
            % Some nessecary data to know.
            %fs = 100; % Sample frequency  Hz.
            
            % Preprocessing of the coordinate file
            
            % Find stationary points
            % p = static points,st = start of movement id, ed = end of movement id
            [p,st,ed] = calibclass.fndzvp(qq);
            % % Check pf what was found
            
            if length(st) ~= length(ed)
                error('Non equal length of start and end pts');
            end
            
            % Region were robot waits
            dWaitreg = min(diff([ed(1:end-1),st(2:end)],1,2))-30;
            
            % if length(st) ~= 150 || length(ed) ~= 150
            %     error('Insufficient nr of start or end points ');
            % end
            if length(st) == 151
                p = p(:,:,2:151);
                ed = ed(2:151);
                %                 st = st(2:151);
            end
            
            % Evaluations for position stabilisation time
            
            % Find entries to P1
            p_ent_ids  = cell(5,1);
            for i = 1:5
                p_ent_ids{i} = ed(6-i:5:end);
            end
            % Center points
            dCenterPoints = mean(p,1);
            dCenterPoints = permute(dCenterPoints,[3 2 1]);
            dCenterPoints = permute(reshape(dCenterPoints',3,5,[]),[2,1,3]);
            
            varargout{1} = [];
            varargout{2} = [];
            % Plot the first point of interest
            for point = 1:5
                fh = figure; set(fh,'name',sprintf('Point nr: %d',point));
                PSt = [];
                PSo = [];
                for j = 1:3
                    %ids = p_ent_ids{point}(j)+[0:499];
                    ids = p_ent_ids{point}(j)+(-30:dWaitreg);
                    subplot(2,2,j);
                    dists = [];
                    for i = ids
                        dists=cat(1,dists,sqrt(sum(diff([qq(i,:);dCenterPoints(6-point,:,j)]).^2)));
                    end
                    plot([1:length(dists)]/calibclass.fs,dists);
                    % This part need to be sorted out, repetability or 0.1 and 0.5. Ask Lars S
                    % theLim = rps(1);
                    %     theLim = 0.1;
                    tmp = dists>theLim;
                    % tmp = dists>0.1;
                    fst = find(gradient(double(tmp))==-0.5,1,'first');
                    rfst = interp1(dists(fst+[0 1]),fst+[0;1],theLim);
                    lst = find(gradient(double(tmp))==-0.5,1,'last')-1;
                    rlst = interp1(dists(lst+[0 1]),lst+[0;1],theLim);
                    if (lst > fst)
                        % Find max deviation from endpoint. "Position overshot" in ISO standard
                        PSo(j) = max(dists(fst:lst));
                        PSt(j) = (rlst-rfst)/calibclass.fs;
                        title(['PSt = ' num2str(PSt(j)) '(s)   PSo = ' num2str(PSo(j)) '(mm)' ])
                    else
                        PSo(j) = 0;
                        PSt(j) = (rlst-rfst)/calibclass.fs;
                        title(['PSt = ' num2str(PSt(j)) '(s)'])
                    end
                    %title(['PSt = ' num2str((lst-fst)/calibclass.fs) '(s)'])
                    hold on
                    ylim = 1.0;
                    % Horisontal line
                    line([1 500]/calibclass.fs,theLim*[1 1],'color','r')
                    line([rfst rlst;rfst rlst]/calibclass.fs,[0 0;ylim ylim],'color','r');
                    hold off
                    axis([0 550/calibclass.fs 0 ylim]);
                    xlabel('Time (s)');
                    ylabel('Dist (mm)');
                end
                
                varargout{1} = cat(2,varargout{1},mean(PSt));
                varargout{2} = cat(2,varargout{2},mean(PSo));
            end
        end
        function varargout=pst(raw,theLim)
            %theLim = 0.1;
            fs  = 100;
            a = [];
            b = [];
            [p,~,ed,offs]=calibclass.fndzvp(raw);
            
            n = length(ed);
            %pp1 = raw(ed(n)+offs,:);
            
            cps = permute(mean(p),[3 2 1]);
            
            jdist = @(q) sqrt(sum(diff(q).^2,2));
            figure;
            for j = 1:n
                dists = jdist([cps(j,:);raw(ed(j)+[0:offs(1)],:)]);
                
                tmp = dists>theLim;
                fst = find(gradient(double(tmp))==-0.5,1,'first');
                rfst = interp1(dists(fst+[0 1]),fst+[0;1],theLim);
                lst = find(gradient(double(tmp))==-0.5,1,'last')-1;
                rlst = interp1(dists(lst+[0 1]),lst+[0;1],theLim);
                if (lst > fst)
                    % Find max deviation from endpoint. "Position overshot" in ISO standard
                    PSo = max(dists(fst:lst));
                    PSt = (rlst-rfst)/fs;
                    titlestr = ['PSt = ' num2str(PSt) '(s)   PSo = ' num2str(PSo) '(mm)' ];
                else
                    PSo = NaN;
                    PSt = (rlst-rfst)/fs;
                    titlestr=['PSt = ' num2str(PSt) '(s)'];
                end
                if true
                    subplot(round(n/5),5,j);
                    plot([1:length(dists)]/fs,dists);
                    hold on
                    ylim = 0.2;
                    % Horisontal line
                    line([1 500]/fs,theLim*[1 1],'color','r')
                    % Two vertical lines
                    line([rfst rlst;rfst rlst]/fs,[0 0;ylim ylim],'color','r');
                    hold off
                    axis([0 550/fs 0 theLim]);
                    xlabel('Time (s)');
                    ylabel('Distance (mm)');
                    title(titlestr);
                end
                a = cat(1,a,PSt);
                b = cat(1,b,PSo);
            end
            varargout{1} = a;
            varargout{2} = b;
            %varargout = {PSt PSo rfst rlst};
        end
        function [maxtab, mintab]=peakdet(v, delta, x)
            %PEAKDET Detect peaks in a vector
            %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
            %        maxima and minima ("peaks") in the vector V.
            %        MAXTAB and MINTAB consists of two columns. Column 1
            %        contains indices in V, and column 2 the found values.
            %
            %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
            %        in MAXTAB and MINTAB are replaced with the corresponding
            %        X-values.
            %
            %        A point is considered a maximum peak if it has the maximal
            %        value, and was preceded (to the left) by a value lower by
            %        DELTA.
            % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
            % This function is released to the public domain; Any use is allowed.
            
            maxtab = [];
            mintab = [];
            
            v = v(:); % Just in case this wasn't a proper vector
            
            if nargin < 3
                x = (1:length(v))';
            else
                x = x(:);
                if length(v)~= length(x)
                    error('Input vectors v and x must have same length');
                end
            end
            
            if (length(delta(:)))>1
                error('Input argument DELTA must be a scalar');
            end
            
            if delta <= 0
                error('Input argument DELTA must be positive');
            end
            
            mn = Inf; mx = -Inf;
            mnpos = NaN; mxpos = NaN;
            
            lookformax = 1;
            
            for i=1:length(v)
                this = v(i);
                if this > mx, mx = this; mxpos = x(i); end
                if this < mn, mn = this; mnpos = x(i); end
                
                if lookformax
                    if this < mx-delta
                        maxtab = [maxtab ; mxpos mx];
                        mn = this; mnpos = x(i);
                        lookformax = 0;
                    end
                else
                    if this > mn+delta
                        mintab = [mintab ; mnpos mn];
                        mx = this; mxpos = x(i);
                        lookformax = 1;
                    end
                end
            end
            
        end
        function varargout = fndzvp(varargin)
            % Find zero velocity points, also displays a boll plot of the detected
            % points.
            % % Get stationary points
            % Indata (in)
            % ********************************************************
            % jaDebug = true;
            
            %fs = 100;               % Sample frequency (Hz)
            st = 40;                % Speed threshold (mm/s)
            staticthresh = false;   % Static threshold flag
            sts = 70;               % Speed threshold (mm/s), static approach
            %cwt = 5;                % Test cycle wait time
            fwt = 3.5;              % Time to wait for robot stabilisation
            nstc = 20; 				% Number of samples to consider
            jaDebug = true;
            refine = false;         % Refine the definition of the startpoint
            
            
            switch nargin
                case 1
                    in = deal(varargin{:});
                case 2
                    % Means extra fine detection of start point, used in positioning
                    % time calculus, set to True if to be used.
                    [in,refine] = deal(varargin{:});
                case 3
                    % true use static threshold for end of movement detection used in
                    % positioning time calculus
                    [in,refine,staticthresh] = deal(varargin{:});
                otherwise
                    error('Wrong number of input arguments');
            end
            
            
            % ********************************************************
            % Special function to find speed threshold
            %
            % jv = @(q,fs) sqrt(sum(diff(q).^2,2))*fs;
            % v = jv(in,100);
            % [a,b] = hist(v);
            % %[c,d]=sort(a);
            % %jojo = sort(b(d(end-1:end)));
            % ste = 0.05*mean(b);
            ste = 10;
            % st = mean(b); % Speed threshold (mm/s)
            
            % ********************************************************
            
            % % Calculation of velocities non filtered
            
            qn = in;
            dS = diff(qn);                  % Moovement
            mvf = sqrt(sum(dS.^2,2))*calibclass.fs;    % Velosity
            % mvf = filter([1 1 1 1 1 1 1],7,mvf); % Filtered velocity
            
            % Find start points
            %ids = mvf<st;
            ids=slusk(in);
            % idss = mvf<sts;
            
            % Can be used to find parts were the robot moves ...
            startpoints = find(gradient(double(ids))==-0.5);
            startpoints = startpoints(1:2:end);
            lbt = 50; % Lenth of back tracking in points
            if refine
                % Refining the the search of the startpoint
                for i = 1:length(startpoints)
                    slask = mvf(startpoints(i)+[-(lbt-1):0]);
                    id  = find(slask < mean(slask(1:5))*1.7,1,'last');
                    startpoints(i) = startpoints(i)-(lbt-id);
                end
            end
            if staticthresh
                endpoints = find(gradient(double(ids))==0.5);
                refine = false;
            else
                endpoints = find(gradient(double(ids))==0.5);
            end
            endpoints = endpoints(1:2:end);
            if refine
                % Refining the establishing of the endpoint
                for i = 1:length(endpoints)
                    id  = find(mvf(endpoints(i)+[0:99]) < ste,1,'first');
                    endpoints(i) = endpoints(i)+id;
                end
            end
            
            if 1
                % Visualizing the selected points
                figure(30);
                t = [1:length(mvf)]/calibclass.fs;
                clf
                plot(t,mvf);grid;xlabel('Time (s)');ylabel('Speed (mm/s)');
                hold on;
                plot(t(endpoints),mvf(endpoints),'ro','MarkerFaceColor','r');
                plot(t(startpoints),mvf(startpoints),'go','MarkerFaceColor','g');
                hold off;
            end
            
            % Number samples jump from calculated stop to stable position
            nBids = round(calibclass.fs*fwt);
            
            % The ids of interest offset
            IdsOffs = nBids+[0:nstc-1]';
            
            % Hitta på ett smartare sätt att bestämma vart mätningen står stilla.
            % Försök nyttja infon om hur långt det är mellan "end" och "start".
            % Bryt ut de sista 85-95% av samplen.
            %
            
            statPoints = [];
            try
                for i = 1:length(endpoints)
                    statPoints = cat(3,statPoints,in(endpoints(i)+IdsOffs,:));
                    if 1
                        hold on
                        plot(t(endpoints(i)+IdsOffs),mvf(endpoints(i)+IdsOffs),'yo','MarkerFaceColor','y');
                        hold off
                    end
                end
            catch
                msgbox('Last point discarded');
                statPoints = [];                % Fault changed 2017-09-25, this row added
                for i = 1:length(endpoints)-1
                    statPoints = cat(3,statPoints,in(endpoints(i)+IdsOffs,:));
                end
            end
            
            if jaDebug
                figure(31);
                mv = permute(mean(statPoints,1),[3,2,1]);
                % 	bollplot(mv(:,1),mv(:,2),mv(:,3));%,repmat([1 0 0],length(mv),1));
                ballplot(mv(:,1),mv(:,2),mv(:,3),'r',20);
                camlight;
                axis vis3d;
            end
            
            varargout{1} = statPoints;
            varargout{2} = startpoints;
            varargout{3} = endpoints;
            varargout{4} = IdsOffs;
            function ut=slusk(qq)
                %qq=getLeicaMeas;
                jaDebug = false;
                % ********* MY MACROS ******************************
                % Macro the calculates the velocity
                % jv = @(q,calibclass.fs) sqrt(sum(diff(q).^2,2))*calibclass.fs;
                %jv = @(q,fs) medfilt1(sqrt(sum(diff(q).^2,2))*calibclass.fs,5);
                %jdist = @(q) sqrt(sum(diff(q).^2,2));
                
                %Plots
                
                % Performs a scatter plot
                %js = @(x,c) scatter3(x(:,1),x(:,2),x(:,3),c);
                %jl = @(x,c) line(x(:,1),x(:,2),x(:,3));
                %jbp = @(x) ballplot(x(:,1),x(:,2),x(:,3),'r',20);
                %cl = @() camlight; % Usage cl()
                % **************************************************
                v = calibclass.jv(qq,100);
                [a,~]=calibclass.peakdet(v,80);
                %
                if jaDebug
                    figure(10)
                    clf
                    subplot(211)
                    plot(a(:,1),a(:,2),'ro')
                    hold on
                    plot(b(:,1),b(:,2),'go')
                    stairs(v)
                end
                % Trial for dynamic thresholding
                npoints = length(v);
                utl=ones(npoints,1); % upper threshold line
                utl(1:a(1,1))=a(1,2);
                for i = 1:length(a)-1
                    utl(a(i,1)+1:a(i+1,1))=linspace(a(i,2),a(i+1,2),diff(a(i+[0 1],1)));
                end
                utl(a(end,1):end)=a(end,2);
                if jaDebug
                    stairs(utl/2,'k')
                    subplot(212)
                    plot(v<utl/2)
                end
                ut = v<utl/2;
            end
        end
        function varargout = fndzvp_v2(varargin)
            % Find zero velocity points, also displays a boll plot of the detected
            % points.
            %% Get stationary points
            % Indata (in)
            % ********************************************************
            % jaDebug = true;
            
            %fs = 100;               % Sample frequency (Hz)
            st = 40;                % Speed threshold (mm/s)
            staticthresh = false;   % Static threshold flag
            sts = 70;               % Speed threshold (mm/s), static approach
            %cwt = 5;                % Test cycle wait time
            fwt = 3.5;              % Time to wait for robot stabilisation
            nstc = 20; 				% Number of samples to consider
            jaDebug = true;
            refine = false;         % Refine the definition of the startpoint
            
            
            switch nargin
                case 1
                    in = deal(varargin{:});
                case 2
                    % Means extra fine detection of start point, used in positioning
                    % time calculus, set to True if to be used.
                    [in,refine] = deal(varargin{:});
                case 3
                    % true use static threshold for end of movement detection used in
                    % positioning time calculus
                    [in,refine,staticthresh] = deal(varargin{:});
                otherwise
                    error('Wrong number of input arguments');
            end
            
            
            % ********************************************************
            % Special function to find speed threshold
            %
            % jv = @(q,calibclass.fs) sqrt(sum(diff(q).^2,2))*calibclass.fs;
            % v = jv(in,100);
            % [a,b] = hist(v);
            % %[c,d]=sort(a);
            % %jojo = sort(b(d(end-1:end)));
            % ste = 0.05*mean(b);
            ste = 10;
            % st = mean(b); % Speed threshold (mm/s)
            
            % ********************************************************
            
            %% Calculation of velocities non filtered
            
            qn = in;
            dS = diff(qn);                  % Moovement
            mvf = sqrt(sum(dS.^2,2))*calibclass.fs;    % Velosity
            % mvf = filter([1 1 1 1 1 1 1],7,mvf); % Filtered velocity
            
            % Find start points
            %ids = mvf<st;
            ids=slusk(in);
            % idss = mvf<sts;
            
            % Can be used to find parts were the robot moves ...
            startpoints = find(gradient(double(ids))==-0.5);
            startpoints = startpoints(1:2:end);
            lbt = 50; % Lenth of back tracking in points
            if refine
                % Refining the the search of the startpoint
                for i = 1:length(startpoints)
                    slask = mvf(startpoints(i)+[-(lbt-1):0]);
                    id  = find(slask < mean(slask(1:5))*1.7,1,'last');
                    startpoints(i) = startpoints(i)-(lbt-id);
                end
            end
            if staticthresh
                endpoints = find(gradient(double(ids))==0.5);
                refine = false;
            else
                endpoints = find(gradient(double(ids))==0.5);
            end
            endpoints = endpoints(1:2:end);
            if refine
                % Refining the establishing of the endpoint
                for i = 1:length(endpoints)
                    id  = find(mvf(endpoints(i)+[0:99]) < ste,1,'first');
                    endpoints(i) = endpoints(i)+id;
                end
            end
            
            if jaDebug
                % Visualizing the selected points
                figure(30);
                t = [1:length(mvf)]/calibclass.fs;
                clf
                plot(t,mvf);grid;xlabel('Time (s)');ylabel('Speed (mm/s)');
                hold on;
                plot(t(endpoints),mvf(endpoints),'ro','MarkerFaceColor','r');
                plot(t(startpoints),mvf(startpoints),'go','MarkerFaceColor','g');
                hold off;
            end
            
            % Number samples jump from calculated stop to stable position
            nBids = round(calibclass.fs*fwt);
            
            % The ids of interest offset
            IdsOffs = nBids+[0:nstc-1]';
            
            % Hitta på ett smartare sätt att bestämma vart mätningen står stilla.
            % Försök nyttja infon om hur långt det är mellan "end" och "start".
            % Bryt ut de sista 85-95% av samplen.
            %
            
            dst = [];
            for i = 2:length(endpoints),
                dst(end+1)=startpoints(i)-endpoints(i-1);
            end
            assert(jrange(dst)/mean(dst)<0.1,'To big time variance between points');
            DST = mean(dst);
            newst = round(0.80*DST);
            newed = round(0.90*DST);
            
            IdsOffs = (newst:newed);
            
            %
            
            statPoints = [];
            try
                for i = 1:length(endpoints)
                    statPoints = cat(3,statPoints,in(endpoints(i)+IdsOffs,:));
                    if jaDebug
                        hold on
                        plot(t(endpoints(i)+IdsOffs),mvf(endpoints(i)+IdsOffs),'yo','MarkerFaceColor','y');
                        hold off
                    end
                end
            catch
                msgbox('Last point discarded');
                statPoints = [];    % Changed 2017-09-25, this row added
                for i = 1:length(endpoints)-1
                    statPoints = cat(3,statPoints,in(endpoints(i)+IdsOffs,:));
                end
            end
            
            if false
                figure(31);
                mv = permute(mean(statPoints,1),[3,2,1]);
                % 	bollplot(mv(:,1),mv(:,2),mv(:,3));%,repmat([1 0 0],length(mv),1));
                ballplot(mv(:,1),mv(:,2),mv(:,3),'r',20);
                camlight;
                axis vis3d;
            end
            
            varargout{1} = statPoints;
            varargout{2} = startpoints;
            varargout{3} = endpoints;
            function ut=slusk(qq)
                %qq=getLeicaMeas;
                %jaDebug = false;
                % ********* MY MACROS ******************************
                % Macro the calculates the velocity
                % jv = @(q,calibclass.fs) sqrt(sum(diff(q).^2,2))*calibclass.fs;
                % jv = @(q,fs) medfilt1(sqrt(sum(diff(q).^2,2))*calibclass.fs,5);
                % jdist = @(q) sqrt(sum(diff(q).^2,2));
                
                %Plots
                
                % Performs a scatter plot
                %js = @(x,c) scatter3(x(:,1),x(:,2),x(:,3),c);
                %jl = @(x,c) line(x(:,1),x(:,2),x(:,3));
                %jbp = @(x) ballplot(x(:,1),x(:,2),x(:,3),'r',20);
                %cl = @() camlight; % Usage cl()
                % **************************************************
                v = calibclass.jv(qq,100);
                
                % ---------- MODIFIED TO SUITE CHINA DATA ---------
                %[a,~]=peakdet(v,80);
                [a,~]=calibclass.peakdet(v,30);
                % ---------- END CHINA MODIFICATION ---------------
                
                %
                if true
                    figure(10)
                    clf
                    subplot(211)
                    plot(a(:,1),a(:,2),'ro')
                    hold on
                    %plot(b(:,1),b(:,2),'go')
                    stairs(v)
                end
                % Trial for dynamic thresholding
                npoints = length(v);
                utl=ones(npoints,1); % upper threshold line
                utl(1:a(1,1))=a(1,2);
                for i = 1:length(a)-1
                    utl(a(i,1)+1:a(i+1,1))=linspace(a(i,2),a(i+1,2),diff(a(i+[0 1],1)));
                end
                utl(a(end,1):end)=a(end,2);
                if false
                    stairs(utl/2,'k')
                    subplot(212)
                    plot(v<utl/2)
                end
                ut = v<utl/2;
            end
            
            function y = jrange(x,dim)
                %JRANGE  Sample range.
                %   Y = JRANGE(X) returns the range of the values in X.  For a vector input,
                %   Y is the difference between the maximum and minimum values.  For a
                %   matrix input, Y is a vector containing the range for each column.  For
                %   N-D arrays, JRANGE operates along the first non-singleton dimension.
                %
                %   JRANGE treats NaNs as missing values, and ignores them.
                %
                %   Y = JRANGE(X,DIM) operates along the dimension DIM.
                %
                %   See also IQR, MAD, MAX, MIN, STD.
                
                
                if nargin < 2
                    y = max(x) - min(x);
                else
                    y = max(x,[],dim) - min(x,[],dim);
                end
            end
        end  
        function h=ballplot(varargin)
            %BALLPLOT 3-D spheres plot.
            %   BALLPLOT(X,Y,Z,C,R) displays colored spheres at the locations specified
            %   by the matrices X,Y,Z (which must all be the same size). The colors of
            %   each sphere are based on the values in C and the radius of each sphere
            %   is determined by the values in R (in axis coordinates).  R can be a
            %   scalar, in which case all the spheres are drawn the same size, or a
            %   vector the same length as prod(size(X)).
            %
            %   When C is a vector the same length as prod(size(X)), the values in C
            %   are linearly mapped to the colors in the current colormap.
            %   When C is a size(X)-by-3 matrix, the values in C specify the
            %   colors of the markers as RGB values.  C can also be a color string.
            %
            %   BALLPLOT(X,Y,Z) draws the spheres with the default size and color.
            %   BALLPLOT(X,Y,Z,C) draws the spheres with a default size.
            %
            %   BALLPLOT(X,Y,Z,C,R,F) draws the spheres having a roundness F. Increase
            %   F in integer steps if the spheres looks crude. Default 1.
            %
            %   BALLPLOT(AX,...) plots into AX instead of GCA.
            %
            %   H = BALLPLOT(...) returns handles to scatter objects created.
            %
            %   Use PLOT3 for single color, single marker size 3-D scatter plots.
            %
            %   Example
            %      [x,y,z] = sphere(8);
            %      X = [x(:)*.5 x(:)*.75 x(:)];
            %      Y = [y(:)*.5 y(:)*.75 y(:)];
            %      Z = [z(:)*.5 z(:)*.75 z(:)];
            %      S = repmat([1 .75 .5]*10,prod(size(x)),1);
            %      C = repmat([1 2 3],prod(size(x)),1);
            %      ballplot(X,Y,Z,C(:));
            %      view(-60,10)
            %      camlight
            %
            %   See also SCATTER3, PLOT3.
            
            [cax,args,nargs] = axescheck(varargin{:});
            narginchk(3,6);
            cax = newplot(cax);
            r=[];
            level=1;
            switch (nargs)
                case 3
                    [x,y,z] = deal(args{1:nargs});
                    [~,c,~] = nextstyle(cax);
                case 4
                    [x,y,z,c] = deal(args{1:nargs});
                case 5
                    [x,y,z,c,r] = deal(args{1:nargs});
                case 6
                    [x,y,z,c,r,level] = deal(args{1:nargs});
                otherwise
                    error('jwtools:ballplot:invalidInput',...
                        'Wrong number of input arguments.');
            end
            % Verify {X,Y,Z}data is correct size
            if any([ndims(y) ndims(z)] ~= ndims(x)) || any([numel(y) numel(z)] ~= numel(x))
                error('jwtools:ballplot:invalidData',...
                    'X,Y, and Z must be vectors of the same size.');
            end
            
            x=x(:);
            y=y(:);
            z=z(:);
            if numel(c)>1
                szc=size(c);
                c=reshape(c,[prod(szc(1:end-1)) szc(end)]);
            end
            % Map colors into colormap colors if necessary.
            if isstr(c)
                [c,msg]=ColorSpecToRGB(c);
                error(msg);
            end
            if isequal(size(c),[1 3]) % string color or scalar rgb
                C = repmat(c,length(x),1);
            elseif length(c)==numel(c) && length(x)==length(c) % is C a vector of size(x)?
                C=colormap(cax);
                mn=min(c);mx=max(c);
                if mx-mn==0
                    C=C(round(length(C)/2),:);
                else
                    C=C(floor((c-mn)/(mx-mn)*(length(C)-1))+1,:);
                end
            elseif isequal(size(c),[length(x) 3]) % vector of rgb's
                C = c;
            else
                error(['C must be a single color, a vector the same length as X, ',...
                    'or a size(X)-by-3 matrix.'])
            end
            
            if isempty(r)
                if ishold
                    ax=axis;
                else
                    hf=figure('Visible', 'off');
                    plot3(x,y,z,'.');
                    axis tight;
                    ax=axis;
                    delete(hf);
                end
                % Get length of shortest axis and make the radius 3% of that
                r=0.03*min(diff(reshape(ax,2,3)));
                r=r*ones(size(x));
            end
            r=r(:);
            sz=size(C);
            C=reshape(C,prod(sz(1:end-1)),3);
            if size(C,1)==1
                C=repmat(C,[size(x,1) 1 1]);
            end
            cax = newplot;
            next = lower(get(cax,'NextPlot'));
            hold_state = ishold(cax);
            np=size(x,1);
            fv=sphere_tri(level, 1);
            
            nv=size(fv.vertices,1);
            nf=size(fv.faces, 1);
            v=repmat(fv.vertices,[1 1 np]);
            if numel(r)==1
                r=repmat(r,[np 1]);
            end
            v=v.*repmat(reshape(r,[1 1 np]), [nv 3 1]);
            v=v+repmat(permute([x y z],[3 2 1]),[nv 1]);
            f=repmat(fv.faces, [1 1 np]) + repmat(reshape(nv*(0:np-1), [1 1 np]), [nf 3 1]);
            c=repmat(permute(C,[3 2 1]), [nf 1]);
            v = reshape(permute(v,[1 3 2]),[],3);
            f = reshape(permute(f,[1 3 2]),[],3);
            c = reshape(permute(c,[1 3 2]),[],3);
            
            hh=patch( ...
                'vertices', v, 'faces', f, 'FaceVertexCData', c ...
                , 'facecolor', 'flat' ...
                , 'edgecolor'	, 'none' ...
                );
            if ~hold_state
                set(gcf,'Renderer', 'OpenGL');
                axis equal tight
                if any(diff(z))
                    view(3)
                    axis vis3d
                else
                    view(2)
                end
                grid
                set(cax,'NextPlot',next);
            end
            if nargout==1
                h=hh;
            end
            
            function [color,msg] = ColorSpecToRGB(s)
                color=[];
                msg = [];
                switch s
                    case 'y'
                        color = [1 1 0];
                    case 'm'
                        color = [1 0 1];
                    case 'c'
                        color = [0 1 1];
                    case 'r'
                        color = [1 0 0];
                    case 'g'
                        color = [0 1 0];
                    case 'b'
                        color = [0 0 1];
                    case 'w'
                        color = [1 1 1];
                    case 'k'
                        color = [0 0 0];
                    otherwise
                        msg = 'unrecognized color string';
                end
            end
            function FV = sphere_tri(maxlevel,r)
                % sphere_tri - generate a triangle mesh approximating a sphere
                %
                % Usage: FV = sphere_tri(shape,Nrecurse,r,winding)
                %
                %   Nrecurse is int >= 0, setting the recursions (default 0)
                %
                %   r is the radius of the sphere (default 1)
                %
                %   FV has fields FV.vertices and FV.faces.  The vertices
                %   are listed in clockwise order in FV.faces, as viewed
                %   from the outside in a RHS coordinate system.
                %
                % The function uses recursive subdivision.  The first
                % approximation is an platonic icosahedron.  Each level of
                % refinement subdivides each triangle face by a factor of 4
                % (see also mesh_refine). At each refinement, the vertices are
                % projected to the sphere surface (see sphere_project).
                %
                % A recursion level of 3 or 4 is a good sphere surface, if
                % gouraud shading is used for rendering.
                %
                % The returned struct can be used in the patch command, eg:
                %
                % % create and plot, vertices: [2562x3] and faces: [5120x3]
                % FV = sphere_tri(',4,1);
                % lighting phong; shading interp; figure;
                % patch('vertices',FV.vertices,'faces',FV.faces,...
                %       'facecolor',[1 0 0],'edgecolor',[.2 .2 .6]);
                % axis off; camlight infinite; camproj('perspective');
                %
                % See also: mesh_refine, sphere_project
                %
                
                % $Revision: 1.15 $ $Date: 2004/05/20 22:28:45 $
                
                % Licence:  GNU GPL, no implied or express warranties
                % Jon Leech (leech @ cs.unc.edu) 3/24/89
                % icosahedral code added by Jim Buddenhagen (jb1556@daditz.sbc.com) 5/93
                % 06/2002, adapted from c to matlab by Darren.Weber_at_radiology.ucsf.edu
                % 05/2004, reorder of the faces for the 'ico' surface so they are indeed
                % clockwise!  Now the surface normals are directed outward.  Also reset the
                % default recursions to zero, so we can get out just the platonic solids.
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                % Twelve vertices of icosahedron on unit sphere
                tau = 0.8506508084; % t=(1+sqrt(5))/2, tau=t/sqrt(1+t^2)
                one = 0.5257311121; % one=1/sqrt(1+t^2) , unit sphere
                
                FV.vertices = [
                    tau,  one,    0
                    -tau,  one,    0
                    -tau, -one,    0
                    tau, -one,    0
                    one,   0 ,  tau
                    one,   0 , -tau
                    -one,   0 , -tau
                    -one,   0 ,  tau
                    0 ,  tau,  one
                    0 , -tau,  one
                    0 , -tau, -one
                    0 ,  tau, -one
                    ];
                % Structure for unit icosahedron
                FV.faces = [  5,  8,  9 ;
                    5, 10,  8 ;
                    6, 12,  7 ;
                    6,  7, 11 ;
                    1,  4,  5 ;
                    1,  6,  4 ;
                    3,  2,  8 ;
                    3,  7,  2 ;
                    9, 12,  1 ;
                    9,  2, 12 ;
                    10,  4, 11 ;
                    10, 11,  3 ;
                    9,  1,  5 ;
                    12,  6,  1 ;
                    5,  4, 10 ;
                    6, 11,  4 ;
                    8,  2,  9 ;
                    7, 12,  2 ;
                    8, 10,  3 ;
                    7,  3, 11 ];
                
                % -----------------
                % refine the starting shapes with subdivisions
                if maxlevel
                    % Subdivide each starting triangle (maxlevel) times
                    for thelevel = 1:maxlevel,
                        
                        % Subdivide each triangle and normalize the new points thus
                        % generated to lie on the surface of a sphere radius r.
                        FV = mesh_refine_tri4(FV);
                        FV.vertices = sphere_project(FV.vertices,r);
                        
                        % An alternative might be to define a min distance
                        % between vertices and recurse or use fminsearch
                        
                    end
                end
            end
            function [ FV ] = mesh_refine_tri4(FV)
                
                % mesh_refine_tri4 - creates 4 triangle from each triangle of a mesh
                %
                % [ FV ] = mesh_refine_tri4( FV )
                %
                % FV.vertices   - mesh vertices (Nx3 matrix)
                % FV.faces      - faces with indices into 3 rows
                %                 of FV.vertices (Mx3 matrix)
                %
                % For each face, 3 new vertices are created at the
                % triangle edge midpoints.  Each face is divided into 4
                % faces and returned in FV.
                %
                %        B
                %       /\
                %      /  \
                %    a/____\b       Construct new triangles
                %    /\    /\       [A,a,c]
                %   /  \  /  \      [a,B,b]
                %  /____\/____\     [c,b,C]
                % A	     c	   C    [a,b,c]
                %
                % It is assumed that the vertices are listed in clockwise order in
                % FV.faces (A,B,C above), as viewed from the outside in a RHS coordinate
                % system.
                %
                % See also: mesh_refine, sphere_tri, sphere_project
                %
                
                
                % ---this method is not implemented, but the idea here remains...
                % This can be done until some minimal distance (D) of the mean
                % distance between vertices of all triangles is achieved.  If
                % no D argument is given, the function refines the mesh once.
                % Alternatively, it could be done until some minimum mean
                % area of faces is achieved.  As is, it just refines once.
                
                
                % $Revision: 1.12 $ $Date: 2004/05/10 21:01:55 $
                
                % Licence:  GNU GPL, no implied or express warranties
                % History:  05/2002, Darren.Weber_at_radiology.ucsf.edu, created
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % NOTE
                % The centroid is located one third of the way from each vertex to
                % the midpoint of the opposite side. Each median divides the triangle
                % into two equal areas; all the medians together divide it into six
                % equal parts, and the lines from the median point to the vertices
                % divide the whole into three equivalent triangles.
                
                % Each input triangle with vertices labelled [A,B,C] as shown
                % below will be turned into four new triangles:
                %
                % Make new midpoints
                % a = (A+B)/2
                % b = (B+C)/2
                % c = (C+A)/2
                %
                %        B
                %       /\
                %      /  \
                %    a/____\b       Construct new triangles
                %    /\    /\       [A,a,c]
                %   /  \  /  \      [a,B,b]
                %  /____\/____\     [c,b,C]
                % A	     c	   C    [a,b,c]
                %
                
                % Initialise a new vertices and faces matrix
                Nvert = size(FV.vertices,1);
                Nface = size(FV.faces,1);
                V2 = zeros(Nface*3,3);
                F2 = zeros(Nface*4,3);
                
                for ff = 1:Nface,
                    
                    % Get the triangle vertex indices
                    NA = FV.faces(ff,1);
                    NB = FV.faces(ff,2);
                    NC = FV.faces(ff,3);
                    
                    % Get the triangle vertex coordinates
                    A = FV.vertices(NA,:);
                    B = FV.vertices(NB,:);
                    CC = FV.vertices(NC,:);
                    
                    % Now find the midpoints between vertices
                    a = (A + B) ./ 2;
                    b = (B + CC) ./ 2;
                    cc = (CC + A) ./ 2;
                    
                    % Find the length of each median
                    %A2blen = sqrt ( sum( (A - b).^2, 2 ) );
                    %B2clen = sqrt ( sum( (B - c).^2, 2 ) );
                    %C2alen = sqrt ( sum( (C - a).^2, 2 ) );
                    
                    % Store the midpoint vertices, while
                    % checking if midpoint vertex already exists
                    [FV, Na] = mesh_find_vertex(FV,a);
                    [FV, Nb] = mesh_find_vertex(FV,b);
                    [FV, Nc] = mesh_find_vertex(FV,cc);
                    
                    % Create new faces with orig vertices plus midpoints
                    F2(ff*4-3,:) = [ NA, Na, Nc ];
                    F2(ff*4-2,:) = [ Na, NB, Nb ];
                    F2(ff*4-1,:) = [ Nc, Nb, NC ];
                    F2(ff*4-0,:) = [ Na, Nb, Nc ];
                    
                end
                
                % Replace the faces matrix
                FV.faces = F2;
                
                return
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [FV, N] = mesh_find_vertex(FV,vertex)
                
                Vn = size(FV.vertices,1);
                Va = repmat(vertex,Vn,1);
                Vexist = find( FV.vertices(:,1) == Va(:,1) & ...
                    FV.vertices(:,2) == Va(:,2) & ...
                    FV.vertices(:,3) == Va(:,3) );
                if Vexist,
                    if size(Vexist) == [1,1],
                        N = Vexist;
                    else,
                        mfvmsg = sprintf('replicated vertices');
                        error(mfvmsg);
                    end
                else
                    FV.vertices(end+1,:) = vertex;
                    N = size(FV.vertices,1);
                end
                
                return
            end
            function V = sphere_project(v,r,c)
                
                % sphere_project - project point X,Y,Z to the surface of sphere radius r
                %
                % V = sphere_project(v,r,c)
                %
                % Cartesian inputs:
                % v is the vertex matrix, Nx3 (XYZ)
                % r is the sphere radius, 1x1 (default 1)
                % c is the sphere centroid, 1x3 (default 0,0,0)
                %
                % XYZ are converted to spherical coordinates and their radius is
                % adjusted according to r, from c toward XYZ (defined with theta,phi)
                %
                % V is returned as Cartesian 3D coordinates
                %
                
                % $Revision: 1.8 $ $Date: 2004/03/29 21:15:36 $
                
                % Licence:  GNU GPL, no implied or express warranties
                % History:  06/2002, Darren.Weber_at_radiology.ucsf.edu, created
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if ~exist('v','var'),
                    spmsg = sprintf('SPHERE_PROJECT: No input vertices (X,Y,Z)\n');
                    error(spmsg);
                end
                
                X = v(:,1);
                Y = v(:,2);
                Z = v(:,3);
                
                if ~exist('c','var'),
                    xo = 0;
                    yo = 0;
                    zo = 0;
                else
                    xo = c(1);
                    yo = c(2);
                    zo = c(3);
                end
                
                if ~exist('r','var'), r = 1; end
                
                % alternate method is to use unit vector of V
                % [ n = 'magnitude(V)'; unitV = V ./ n; ]
                % to change the radius, multiply the unitV
                % by the radius required.  This avoids the
                % use of arctan functions, which have branches.
                
                
                % Convert Cartesian X,Y,Z to spherical (radians)
                theta = atan2( (Y-yo), (X-xo) );
                phi   = atan2( sqrt( (X-xo).^2 + (Y-yo).^2 ), (Z-zo) );
                % do not calc: r = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2);
                
                %   Recalculate X,Y,Z for constant r, given theta & phi.
                R = ones(size(phi)) * r;
                xx = R .* sin(phi) .* cos(theta);
                yy = R .* sin(phi) .* sin(theta);
                zz = R .* cos(phi);
                
                V = [xx yy zz];
                
                return
            end
            function [ax,args,nargs] = axescheck(varargin)
                %AXESCHECK Process leading Axes object from input list
                %   [AX,ARGS,NARGS] = AXESCHECK(ARG1,ARG2,...) checks if ARG1 is an Axes
                %   and returns it in AX if it is and returns the processed argument
                %   list in ARGS and NARGS.  If ARG1 is not an Axes, AX will return empty.
                %   Also checks arguments that are property-value pairs 'parent',ARG.
                
                %    Copyright 1984-2003 The MathWorks, Inc.
                %    $Revision $  $Date: 2004/04/10 23:28:36 $
                
                args = varargin;
                nargs = nargin;
                ax=[];
                if (nargs > 0) && (numel(args{1}) == 1) && ...
                        ishandle(args{1}) && strcmp(get(args{1},'type'),'axes')
                    ax = args{1};
                    args = args(2:end);
                    nargs = nargs-1;
                end
                if nargs > 0
                    inds = find(strcmpi('parent',args));
                    if ~isempty(inds)
                        inds = unique([inds inds+1]);
                        pind = inds(end);
                        if nargs <= pind && ishandle(args{pind})
                            ax = args{pind};
                            args(inds) = [];
                            nargs = length(args);
                        end
                    end
                end
            end
        end
        function pl3d(in)
            rm = @(x) x - repmat(mean(x),size(x,1),1); % remove mean
            jl = @(x) plot3(x(:,1),x(:,2),x(:,3)); % plot 3d
            fh = findobj('Name','pl3d');
            if isempty(fh)
                fh  = figure;
                fh.Name = 'pl3d';
                hold on;
            else
                figure(fh);
            end
            %             jl(rm(in));
            jl(in);
            ca = gca;
            childs = ca.Children;
            co = ca.ColorOrder;
            for i = 1:length(childs)
                childs(i).Color = co(i,:);
            end
            xlabel('X');ylabel('Y');zlabel('Z');
            grid on;
        end
        function trace(varargin)
            % Function follows the Leica path with a small Dodecahedron. The
            % usage is as follows:
            %   a = calibclass;
            %   raw = a.jload;
            %   a.pl3d(raw)
            %   a.trace
            
            % Parse input
            if nargin == 0
                iterStep = 10;
            elseif nargin == 1
                iterStep = varargin{1};
            else
                error('Not sufficient nr of input arguments..')
            end
            % Create the Dodecahedron
            r=1;
            phi=(1+sqrt(5))/2;
            V1=[1;(1/phi);-phi;phi;-1;0;-phi;1;-1;-1;1;(1/phi);-1;0;0;-(1/phi);phi;-(1/phi);1;0;];
            V2=[1;0;-(1/phi);(1/phi);1;-phi;(1/phi);-1;1;-1;-1;0;-1;-phi;phi;0;-(1/phi);0;1;phi;];
            V3=[[1;phi;0;0;-1;-(1/phi);0;1;1;1;-1;-phi;-1;(1/phi);-(1/phi);phi;0;-phi;-1;(1/phi);]];
            F=[1,2,16,9,20;2,16,10,14,8;16,9,7,3,10;7 9 20 15 5;5,7,3,13,18;3,13,6,14,10;6,13,18,12,11;6,11,17,8,14;11,12,19,4,17;1,2,8,17,4;1,4,19,15,20;12,18,5,15,19];
            [THETA,PHI,R]=cart2sph(V1,V2,V3);
            R=r.*ones(size(V1(:,1)));
            [V1,V2,V3]=sph2cart(THETA,PHI,R);
            V=[V1 V2 V3];
            % Figure out the properties of the current window
            axis equal; grid on; hold on; view(3); %axis off;
            ax = axis;
            range = @(a) abs(max(a)-min(a));
            mxr = max([range(ax(1:2)),range(ax(3:4)),range(ax(5:6))]);
            kk = mxr/40;
            axis(ax + [-1 1 -1 1 -1 1]*kk*2)
            lobj = findobj('type','line');
            xyz = [get(lobj,'xdata')',get(lobj,'ydata')',get(lobj,'zdata')'];
            
            jfixV = @(v,x,s) (v*s)+repmat(x,[size(v,1),1]);
            Vny = jfixV(V,xyz(1,:),kk);
            delete(findobj('Tag','Dodecahedron'))
            ph = patch('Faces',F,'Vertices',Vny,'FaceColor','b','FaceAlpha', ...
                0.6,'EdgeColor','k','LineWidth',2,'Tag','Dodecahedron');
            for i = 2:iterStep:length(xyz)
                Vny = jfixV(V,xyz(i,:),kk);
                set(ph,'Vertices',Vny);
                pause(0.001);
                drawnow;
            end
        end
        function [ut,fn] = jload(varargin)
            if nargin == 2
                fn = calibclass.jgetfile(varargin{:});
            elseif nargin == 1
                fn = varargin{1};
            elseif nargin == 0
                fn = calibclass.jgetfile();
            end
            if fn ~= 0
                fid = fopen(fn,'r');
                in = fread(fid,Inf,'double');
                fclose(fid);
                nsamp = length(in);
                nCords = nsamp/3;
                ut=reshape(in,[3,nCords])';
            end
        end
        function frontview(in)
        end        
        function ut = meanfilt(in,l)
            windowSize = l;
            b = (1/windowSize)*ones(1,windowSize);
            a = 1;
            ut = [];
            for i = 1:3
                ut = cat(2,ut,filter(b,a,in(:,i)));
            end
        end
        function out = rearange(in)
            pm = mean(in);
            out = permute(pm,[3 2 1]);
        end
        function [ap,rp] = aprp(in)
            % APRP  Function calculates point repeatability according to
            % ISO 9283
            %
            % Usage : [ap,rp] = aprp(X)
            % The input argument X requires a n x 3 [x,y,z] matrix of same
            % command point of repetitive measurements to evaluate.
            %
            
            % Calculate centre points
            cp = mean(in);
            % Calculate distances from center point (dfcp)
            dfcp = @(p,pc) sqrt((p(:,1)-pc(1)).^2 + (p(:,2)-pc(2)).^2 + (p(:,3)-pc(3)).^2);
            dp = dfcp(in,cp);
            % RP according to ISO
            frp = @(d) mean(d) + 3*std(d);
            rp = frp(dp);
            % AP according to ISO
            ap = dfcp(in(1,:),mean(in(2:end,:)));
        end
        function varargout = pathaccuracy(qq,varargin)
            % ATRT calculus
            % Usage : [at,rt] = a.pathaccuracy(rawdata);
            %         [at,rt] = a.pathaccuracy(rawdata,tail);
            % raedata : [x,y,z] data in a column array were fs = 100Hz
            % tail : mm distance removed from ends of path. Default 50 mm
            
            % Variables/constants
            if nargin > 1
                tails = varargin{1};
            else
                tails = 50;         % nof mm to cut from start and stop
            end
            jaDebug = true;
            % Find path accurancy for ISO - xxxx
            
            % P4 is the start point, and the path to evaluate is P4 -> P2 ten times
            
            [p,st,ed] = calibclass.fndzvp(qq,true);
            
            if length(st) ~= length(ed)
                error('The number of start and endpoints differ ');
            end
            
            if length(st) < 20
                error('The number of detected points is less than 20 ');
            end
            
            % Remove the heating up at the beginning
            p = p(:,:,end-19:end);
            st = st(end-19:end);
            ed = ed(end-19:end);
            
            
            % Steps trying to evaluate path acurancy
            
            % Formula for the distance from a point to a line, P1 and P2 points
            % defining the line. P0 point to which distance is evaluated.
            js = @(x,c) scatter3(x(:,1),x(:,2),x(:,3),c);
            pdist = @(P0,P1,P2) norm(cross(diff([P1;P2]),diff([P1;P0])))/norm(diff([P1;P2]));
            
            % Trying to find the line segments to evaluate.
            % For this example try P2 to P1
            
            % old ,pths = [starts(5:5:end),endpts(5:5:end)];
            varargout{1} = [];
            varargout{2} = [];
            for direct = 1:2
                switch direct
                    case 1
                        % Sum up the P4 -> P2 movements
                        %             pths = [st(1:2:end), ed(1:2:end)];
                        pths = [st(1:2:end), st(2:2:end)];
                        p4 = cat(4,p(:,:,1:2:end),p(:,:,2:2:end));
                    case 2
                        % Sum up the P2 -> P4 movements
                        %             pths = [st(2:2:end), ed(2:2:end)];
                        tmpSt = st(3:2:end);
                        pths = [st(2:2:end), [tmpSt;tmpSt(end)+median(diff(tmpSt))]];
                        p4 = cat(4,p(:,:,2:2:end),p(:,:,1:2:end));
                end
                
                % Define were to start i.e. 10mm after start
                fds = @(x) sqrt(sum(diff(x,1,1).^2,2));
                
                % Collect all paths
                xyzfl = {};
                for i = 1:10
                    xyzfl{i} = qq(pths(i,1):pths(i,2),:);
                end
                xyz = xyzfl;
                
                if false
                    figure(1)
                    clf(1)
                    %jp(qq)
                    hold on
                    for i = 1:10
                        js(xyz{i},'r')
                    end
                    hold off
                end
                
                
                
                % Define planes for evaluation
                cps = mean(mean(p4,1),3); % Get average centerpoints
                cps = permute(cps,[4 2 1 3]);
                
                % Print the paths in a figure
                if true
                    calibclass.jvr(xyzfl,cps(1,:),cps(2,:),tails);
                    title(sprintf('Direction : %d',direct));
                end
                
                % Calculate the normal vector for line P2 -> P1
                nP2P1 = diff(cps([1 2],:));
                enP2P1 = nP2P1/norm(nP2P1);
                
                % Calculate points on virtual line
                % t = 0:1/(nPointsMin-1):1;
                % vL = [t(:) ones(length(t),1)] * [nP2P1-enP2P1*2*tails ; cps(1,:)+enP2P1.*tails];
                % vL = vL(2:end-1,:);
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % An alternate formulation that divides the distance between the two center
                % points in for instance 10 pieces.
                
                % tails 10mm gives
                cpsnew = [];
                cpsnew(1,:) = cps(1,:)+enP2P1*tails;
                cpsnew(2,:) = cps(2,:)-enP2P1*tails;
                
                % divides the line in 20 segments
                divisor = 20;
                
                x = linspace(cpsnew(1,1),cpsnew(2,1),divisor);
                y = linspace(cpsnew(1,2),cpsnew(2,2),divisor);
                z = linspace(cpsnew(1,3),cpsnew(2,3),divisor);
                vL = [x(:) y(:) z(:)];
                
                if jaDebug
                    %jbp = @(x,c,r) ballplot(x(:,1),x(:,2),x(:,3),c,r);
                    figure(3);clf;
                    bh(1) = calibclass.jbpcr(cps,'r',20);
                    camlight
                    axis vis3d
                    hold on
                    bh(2) = calibclass.jbpcr(vL,'b',10);
                    set(bh(1),'facealpha',0.3);
                    hold off
                end
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %Plane coordinate points, remove line endpoints
                
                
                nofLines = 10;
                nofPlanes = length(vL);
                isps = zeros(nofLines,4,nofPlanes); % Intersection points
                for j = 1:nofPlanes  % Plane number
                    for LNo = 1:nofLines % Line number
                        pos=[];
                        for i = 1:length(xyz{LNo})
                            pos = cat(1,pos,dot(enP2P1,diff([xyz{LNo}(i,:);vL(j,:)])));
                        end
                        mm = pos>0;
                        if mm(1)==1
                            mm = 1 - mm;
                        end
                        pp = find(mm==1,1,'first');
                        pp = [pp-1,pp];
                        if length(pp)>2
                            disp('Error');
                        end
                        [I,check]=calibclass.plane_line_intersect(enP2P1,vL(j,:),xyz{LNo}(pp(1),:),xyz{LNo}(pp(2),:));
                        isps(LNo,:,j) = [I,check];
                    end
                end
                
                if false
                    % Visualisation
                    clf
                    hold on
                    for i = 1:nofPlanes, js(isps(:,:,i),'b');end
                    hold off
                end
                % Mean distances on the planes, and path accurancy
                ATp = [];
                RTp = [];
                for i = 1:nofPlanes
                    ATp = cat(1,ATp,sqrt(sum(diff([vL(i,:);mean(isps(:,1:3,i),1)]).^2)));
                    RTpt =[];
                    for j = 1:length(isps(:,1:3,i))
                        RTpt = cat(1,RTpt,sqrt(sum(diff([isps(j,1:3,i);mean(isps(:,1:3,i),1)]).^2)));
                    end
                    RTp = cat(1,RTp,mean(RTpt)+3*std(RTpt));
                end
                varargout{1} = cat(2,varargout{1},max(ATp));
                varargout{2} = cat(2,varargout{2},max(RTp));
            end
        end
        function varargout = pathveleval(qq,commandVel)
            % PATHEVELVAL  Function evaluates the velocity performance
            %              According to ISO 9283.
            %
            % Usage :  [av,rv,fv]=pathveleval(raw,cV)
            %          raw : [x,y,z] data from Leica logg (mm) fs = 100 Hz
            %          cV  : Comamnd velocity (mm/s)
            
            
            %commandVel = 1600;             % Command velocity (mm/s)
            fs = calibclass.fs;             % Sample frequency (Hz)
            
            % Find path accurancy for ISO - xxxx
            
            % P4 is the start point, and the path to evaluate is P4 -> P2 ten times
            [p,st,ed] = calibclass.fndzvp(qq);
            
            if length(st) ~= length(ed)
                error('The number of start and endpoints differ ');
            end
            
            if length(st) < 20
                error('The number of detected points is less than 20 ');
            end
            
            % Remove the heating up at the beginning
            p = p(:,:,end-19:end);
            st = st(end-19:end);
            ed = ed(end-19:end);
            jp1 = mean(mean(p(:,:,1:2:end),3),1);
            jp2 = mean(mean(p(:,:,2:2:end),3),1);
            lencubediag = sqrt(sum(diff([jp2;jp1]).^2));% - sqrt(3)*1000;
            
            
            % Steps trying to evaluate path acurancy
            
            % Formula for the distance from a point to a line, P1 and P2 points
            % defining the line. P0 point to which distance is evaluated.
            % js = @(x,c) scatter3(x(:,1),x(:,2),x(:,3),c);
            % pdist = @(P0,P1,P2) norm(cross(diff([P1;P2]),diff([P1;P0])))/norm(diff([P1;P2]));
            pdist2 = @(P0,P1) sqrt(sum(diff([P0;P1]).^2));
            % Trying to find the line segments to evaluate.
            % For this example try P2 to P1
            
            % old ,pths = [starts(5:5:end),endpts(5:5:end)];
            varargout{1} = [];
            varargout{2} = [];
            varargout{3} = [];
            % **************************************************************************
            for direct = 1:2
                switch direct
                    case 1
                        % Sum up the P4 -> P2 movements
                        stFigName = 'P4 -> P2';
                        pths = [st(1:2:end), ed(1:2:end)];
                        %p4 = cat(4,p(:,:,1:2:end),p(:,:,2:2:end));
                    case 2
                        % Sum up the P2 -> P4 movements
                        stFigName = 'P2 -> P4';
                        pths = [st(2:2:end), ed(2:2:end)];
                        %p4 = cat(4,p(:,:,2:2:end),p(:,:,1:2:end));
                end
                % Path velocity code here
                % Calculation of velocities non filtered
                %     dqq = diff(qq);
                %     mvnf = sqrt(sum(dqq.^2,2))*fs;
                %jm;
                mvnf = calibclass.jv(qq,fs);
                vxyz = {};
                pathMeanVel = [];
                mimx = [];
                fh = figure(direct+3);
                set(fh,'name',stFigName)
                clf
                hold on
                for i = 1:10
                    subplot(3,4,i);
                    qqq = qq(pths(i,1):pths(i,2),:);
                    vxyz{i} = mvnf(pths(i,1):pths(i,2),:);
                    jst = find(gradient(vxyz{i})<0,1,'first'); % Finds first non accelerating point
                    jed = find(gradient(vxyz{i})>0,1,'last');  % Finds last non accelerating point
                    % Check duration of stable velocity
                    if pdist2(qqq(jst,:),qqq(jed,:))/lencubediag >= 0.5
                        % Duration enough
                        stairs(vxyz{i});%title(['MeanVelocity = ' num2str(mean(vxyz{i})) '(mm/s)']);
                        line([jst jst],[min(vxyz{i}) max(vxyz{i})],'color','r');
                        line([jed jed],[min(vxyz{i}) max(vxyz{i})],'color','r');
                        mimx = cat(1,mimx,[min(vxyz{i}(jst:jed)) max(vxyz{i}(jst:jed))]);
                        pathMeanVel = cat(1,pathMeanVel,mean(vxyz{i}(jst:jed)));
                        title(['MeanVelocity = ' num2str(pathMeanVel(end)) '(mm/s)']);
                    else
                        % To short duration
                        mimx = cat(1,mimx,[NaN NaN]);
                        pathMeanVel = cat(1,pathMeanVel,NaN);
                        title(['MeanVelocity = ' 'NaN' '(mm/s)']);
                    end
                end
                hold off
                
                AV = (mean(pathMeanVel) - commandVel)/commandVel*100;
                RV = 3*std(pathMeanVel)/commandVel*100;
                FV = max(diff(mimx,1,2));
                varargout{1} = cat(2,varargout{1},AV);
                varargout{2} = cat(2,varargout{2},RV);
                varargout{3} = cat(2,varargout{3},FV);
            end
        end
        function [] = driftcheck(points,fignum)
            %DRIFTCHECK Checks APRP raw data for drift
            %   This function checks for drift of measurement in its raw data
            %   displays the result in a coulored scatter plot
            %   Usage : 
            %           a = calibclass
            %           raw = a.jload;         % Load a APRP file 30 laps
            %           a.calibrate            % Optional
            %           rawrc = a.lc2rc(raw);  % Optional
            %           p = a.fndzvp_v2(rawrc);
            %           a.driftcheck(p,1)
            
            if nargin == 0
            else
                s3 = size(points,3);
                % Quality check of inputdata
                if size(points,2)==3 && mod(s3,5)==0
                    % Proper size
                else
                    % Inproper size
                    error('Illegal size of input data, must be [~,3,150]');
                end
            end
            % group in points
            
            id{5} = 1:5:s3;
            id{4} = 2:5:s3;
            id{3} = 3:5:s3;
            id{2} = 4:5:s3;
            id{1} = 5:5:s3;
            p = {};
            CatDir = 3;
            
            for i = 1:5
                tmp = [];
                for ids = id{i}
                    tmp = cat(CatDir,tmp,points(:,:,ids));
                end
                p{i} = tmp;
            end
            
            ColOrd = jet;
            js = @(x,c) scatter3(x(:,1,1),x(:,2,1),x(:,3,1),3,repmat(c,[length(x),1]));
            
            figure(fignum);clf(fignum);
            ax = [];
            h = [];
            for i = 1:5
                subplot(2,3,i)
                hold on;
                ci = linspace(1,64,s3/5);
                for j = 1:s3/5
                    h(i) = js(p{i}(:,:,j),ColOrd(ci(j),:));
                    set(h(i),'Tag',num2str(j));
                end
                axis equal
                grid;
                xlabel('X');
                ylabel('Y');
                zlabel('Z');
                title(['Point : ' num2str(i)]);
                view(3);
                ax = cat(1,ax,axis);
            end
            % Equalise axes
            X = ax(:,1:2);Y = ax(:,3:4);Z = ax(:,5:6);
            dX = max(diff(X,1,2));dY = max(diff(Y,1,2));dZ = max(diff(Z,1,2));
            mX = mean(X,2);mY = mean(Y,2);mZ = mean(Z,2);
            newX = [mX-dX/2,mX+dX/2];newY = [mY-dY/2,mY+dY/2];newZ = [mZ-dZ/2,mZ+dZ/2];
            newAx = [newX,newY,newZ];
            for i = 1:5
                subplot(2,3,i)
                axis(newAx(i,:));
            end
        end
        function dispRobot(sf,pos)
            load robot;
            patch(-sf*y,sf*x,sf*z,c);
            axis equal;
            light;
            shading interp;
        end
    end
end