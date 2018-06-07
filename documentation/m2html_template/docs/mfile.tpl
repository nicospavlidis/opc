<!DOCTYPE html>
<html>
<head>
  <title>Description of {NAME}</title>
  <meta name="keywords" content="{NAME}">
  <meta name="description" content="{H1LINE}">
  <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
  <meta name="generator" content="m2html v1.5 &copy; 2003-2005 Guillaume Flandin">
  <meta name="robots" content="index, follow">
  <link type="text/css" rel="stylesheet" href="{MASTERPATH}m2html.css">
  <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
  <link type="text/css" rel="stylesheet" href="{MASTERPATH}docs.css">
</head>
<body>
<a name="_top"></a>
<!-- BEGIN pathline -->
<!-- -->
<!-- END pathline -->
<div class="wrapper">
	<nav id="sidebar">
		<div class="sidebar-header">
			<h3>Matlab Directories</h3>
        </div>
        <ul class="list-unstyled components">
			<li class="">
				<a href="{MASTERPATH}{INDEX}">Home</a>
            </li>
			<!-- BEGIN rowdir --> 
			<li class="submenu">
				<a href="#{DIR_NAME}" data-toggle="collapse" aria-expanded="false">{DIR_NAME}</a>
                <ul class="collapse list-unstyled" id="{DIR_NAME}">
				<!-- BEGIN row-m -->
					<li class="{NAME}"><a href="{MASTERPATH}{FILE_SOURCE}" title="{H1LINE}">{NAME} <!-- BEGIN mexfile --> <img src="{MASTERPATH}mex.png" alt="MEX" border="0"> <!-- END mexfile --> </a></li>
				<!-- END row-m -->	
                </ul>
            </li>
			<!-- END rowdir -->
        </ul>
    </nav>

	<div id="content">
        <button type="button" id="sidebarCollapse" class="btn btn-info navbar-btn">
			<i class="glyphicon glyphicon-align-left"></i>
            <span>Toggle Menu</span>
        </button>
		<h1>Matlab Index</h1>
		<h1>{NAME}
<!-- BEGIN mexfile --> &nbsp;&nbsp;<img src="{MASTERPATH}{MEXTYPE}.png" alt="{PLATFORMS}" border="0" title="{PLATFORMS}"> <!-- END mexfile -->
</h1>

<h2><a name="_name"></a>PURPOSE <a class="top" href="#_top"><span class="glyphicon glyphicon-arrow-up" aria-hidden="true"></span></a></h2>
<div class="box"><strong>{H1LINE}</strong></div>

<h2><a name="_synopsis"></a>SYNOPSIS <a class="top" href="#_top"><span class="glyphicon glyphicon-arrow-up" aria-hidden="true"></span></a></h2>
<div class="box"><strong>{SYNOPSIS} <!-- BEGIN script --> This is a script file. <!-- END script --> </strong></div>

<h2><a name="_description"></a>DESCRIPTION <a class="top" href="#_top"><span class="glyphicon glyphicon-arrow-up" aria-hidden="true"></span></a></h2>
<div class="fragment"><pre class="comment">{DESCRIPTION}</pre></div>

<!-- crossreference -->
<h2><a name="_cross"></a>CROSS-REFERENCE INFORMATION <a class="top" href="#_top"><span class="glyphicon glyphicon-arrow-up" aria-hidden="true"></span></a></h2>
This function calls:
<ul style="list-style-image:url({MASTERPATH}matlabicon.gif)">
<!-- BEGIN crossrefcall -->
<li><a href="{L_NAME_CALL}" class="code" title="{SYNOP_CALL}">{NAME_CALL}</a>	{H1LINE_CALL}</li>
<!-- END crossrefcall -->
</ul>
This function is called by:
<ul style="list-style-image:url({MASTERPATH}matlabicon.gif)">
<!-- BEGIN crossrefcalled -->
<li><a href="{L_NAME_CALLED}" class="code" title="{SYNOP_CALLED}">{NAME_CALLED}</a>	{H1LINE_CALLED}</li>
<!-- END crossrefcalled -->
</ul>
<!-- crossreference -->

<!-- BEGIN subfunction -->
<h2><a name="_subfunctions"></a>SUBFUNCTIONS <a class="top" href="#_top"><span class="glyphicon glyphicon-arrow-up" aria-hidden="true"></span></a></h2>
<ul style="list-style-image:url({MASTERPATH}matlabicon.gif)">
<!-- BEGIN onesubfunction -->
<li><a href="{L_SUB}" class="code">{SUB}</a></li>
<!-- END onesubfunction -->
</ul>
<!-- END subfunction -->

<!-- BEGIN download -->
<h2><a name="_download"></a>DOWNLOAD <a class="top" href="#_top"><span class="glyphicon glyphicon-arrow-up" aria-hidden="true"></span></a></h2>
<p><a href="{NAME}.m">{NAME}.m</a></p>
<!-- END download -->

<!-- BEGIN source -->
<h2><a name="_source"></a>SOURCE CODE <a class="top" href="#_top"><span class="glyphicon glyphicon-arrow-up" aria-hidden="true"></span></a></h2>
<div class="fragment"><pre>{SOURCECODE}</pre></div>
<!-- END source -->
<hr><address>Generated on {DATE} by <strong><a href="http://www.artefact.tk/software/matlab/m2html/" title="Matlab Documentation in HTML">m2html</a></strong> &copy; 2005</address>
</div>
</div>

<script src="https://code.jquery.com/jquery-1.12.0.min.js"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
<script>
$(document).ready(function () {
	$('.{NAME}').parents('li').addClass('active').contents('ul').addClass('in').attr('aria-expanded',true);
	$('.{NAME} a').css('background', '#5bc0de');
    $('#sidebarCollapse').on('click', function () {
        $('#sidebar, #content').toggleClass('active');
        $('.collapse.in').toggleClass('in');
        $('a[aria-expanded=true]').attr('aria-expanded', 'false');
    });
});
</script>
</body>
</html>
