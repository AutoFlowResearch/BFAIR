import React, { Component } from "react";
import "../../../App.css";

class FooterArea extends Component {
  state = {};
  render() {
    return (
      <footer className="page-footer font-small blue">
        <div className="footer-copyright text-center py-3">
          © 2018 Copyright:
          <a href="https://www.zymr.com/"> MDBootstrap.com</a>
        </div>
      </footer>
    );
  }
}

export default FooterArea;
