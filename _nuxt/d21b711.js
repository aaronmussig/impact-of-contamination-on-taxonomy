(window.webpackJsonp=window.webpackJsonp||[]).push([[1],{200:function(t,r,e){var content=e(276);content.__esModule&&(content=content.default),"string"==typeof content&&(content=[[t.i,content,""]]),content.locals&&(t.exports=content.locals);(0,e(61).default)("4d6a2e0f",content,!0,{sourceMap:!1})},213:function(t,r,e){"use strict";var n=e(305),o=e(307),c=e(306),l={name:"DefaultLayout",data:function(){return{}}},f=e(75),component=Object(f.a)(l,(function(){var t=this._self._c;return t(n.a,[t(c.a,[t(o.a,[t("Nuxt")],1)],1)],1)}),[],!1,null,null,null);r.a=component.exports},221:function(t,r,e){e(222),t.exports=e(223)},275:function(t,r,e){"use strict";e(200)},276:function(t,r,e){var n=e(60)(!1);n.push([t.i,"h1[data-v-35e10596]{font-size:20px}",""]),t.exports=n},57:function(t,r,e){"use strict";var n=e(305),o={name:"EmptyLayout",layout:"empty",props:{error:{type:Object,default:null}},data:function(){return{pageNotFound:"404 Not Found",otherError:"An error occurred"}},head:function(){return{title:404===this.error.statusCode?this.pageNotFound:this.otherError}}},c=(e(275),e(75)),component=Object(c.a)(o,(function(){var t=this,r=t._self._c;return r(n.a,{attrs:{dark:""}},[404===t.error.statusCode?r("h1",[t._v("\n    "+t._s(t.pageNotFound)+"\n  ")]):r("h1",[t._v("\n    "+t._s(t.otherError)+"\n  ")]),t._v(" "),r("NuxtLink",{attrs:{to:"/"}},[t._v("\n    Home page\n  ")])],1)}),[],!1,null,"35e10596",null);r.a=component.exports}},[[221,8,2,9]]]);