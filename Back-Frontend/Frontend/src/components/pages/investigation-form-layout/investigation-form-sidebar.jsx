import React, { Component } from 'react';

export class GeneralinvestigationSidebar extends Component {
  render() {
    return (
      <div id="main">
        <div class="container">
          <div class="accordion" id="faq">
            <div class="card">
              <div class="card-header" id="faqhead1">
                <a
                  href="#"
                  class="btn btn-header-link"
                  data-toggle="collapse"
                  data-target="#sideMenu"
                  aria-expanded="true"
                  aria-controls="faq1"
                >
                  <b>Studdy 1</b>
                </a>
              </div>

              <div
                id="sideMenu"
                class="collapse show"
                aria-labelledby="faqhead1"
                data-parent="#faq"
              >
                <div class="card-body">
                  <p>Assay 1</p>
                  <p>Assay 2</p>
                  <p>Assay 3</p>
                  <p>Assay 4</p>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    );
  }
}

export default GeneralinvestigationSidebar;
