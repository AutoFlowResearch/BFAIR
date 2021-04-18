import React, { Component, Fragment } from "react";
import HeaderArea from "../common/headerArea";
import SidebarArea from "../common/sidebarArea";
import ContentArea from "../common/contentArea";
import FooterArea from "../common/footerArea";
import Container from "../common/Container";
import Row from "../common/Row";
import Column from "../common/Column";

class MainPage extends Component {
  render() {
    return (
      <Fragment>
        <HeaderArea />
        <Container>
          <Row>
            <SidebarArea {...this.props} />
            <ContentArea {...this.props} />
          </Row>
        </Container>
        <div className="control-sidebar-bg" />

        {/* <FooterArea /> */}
      </Fragment>
    );
  }
}

export default MainPage;
